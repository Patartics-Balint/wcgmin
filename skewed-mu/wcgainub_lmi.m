function [gamma, gamw, Dz, Dv] = wcgainub_lmi(G, blk, bases, gamopt, use_kyp)
% function out = wcgainub_lmi(G,blk,bases)
%
% Compute upper bound on worst-case gain using linear combination of
% basis functions for scaling.
% LFT notation is  
%   [z;e] = G [v;d]  
%    v = Delta z
%
% Inputs
%   G := (Nz+Ne)-by-(Nv+Nd) FRD object 
%   blk := Nblk-by-2 matrix describing uncertainty blocks. All
%     uncertainties are assumed to be full complex blocks.
%   bases : = Nblk-by-1 cell array where bases{i} contains a Ni-by-1
%      vector of state-spaces systems used as bases vectors for the
%      scaling associated with the i^th block.
%   [Possibly include optional options input?]
%
% Outputs
%   gamopt := ???
%   D := ???
%

%% Get Dimensions
Nblk = size(blk,1);
Nz = sum(blk(:,2));
Nv = sum(blk(:,1));
Ne = size(G,1) - Nz;
Nd = size(G,2) - Nv;

%% Get Frequency Vector
w = G.Frequency;
Nw = numel(w);

%% Build Bases as SS Objects
Psiv = [];
Psiz = [];
for i=1:Nblk
    Psii = frd(bases{i},w);
    for j = 1:blk(i,1)
        Psiv = blkdiag( Psiv, Psii );
    end    
    %Psiv = blkdiag( Psiv, kron( eye(blk(i,1)), Psii) );
    
    for j = 1:blk(i,2)
        Psiz = blkdiag( Psiz, Psii );
    end    
    %Psiz = blkdiag( Psiz, kron( eye(blk(i,2)), Psii) );
end    


%% Initialize LMI and Create Variables

% Initialize LMI
setlmis([]);

% Create variables for gain at each frequency and peak gain
%   Each gamw(i) represents the scalar gain at frequency i.
%   gamwI(i) represents gamw(i)*I_Ne
%   gamwZI(i) represents blkdiag( zeros(Nw), gam(i)*I_Nd )
%   gamwDiag represents diag( gamw )
gamw = zeros(Nw,1);
gamwI = zeros(Nw,1);
gamwZI = zeros(Nw,1);
for i=1:Nw
    [gamw(i),~,Xdec] = lmivar(1,[1 1]);
    [gamwI(i),~,Xdec] = lmivar(3,i*eye(Ne));
    [gamwZI(i),~,Xdec] = lmivar(3,blkdiag(zeros(Nv),i*eye(Nd)));    
end
gamwDiag = lmivar(3, diag(1:Nw) );

% Create matrix variables gampeak*I
%gampeak = lmivar(1,[1 1]);
gampeakDiag = gamopt * (1 + 1e-1) * eye(Nw);

% Create variables for scaling
lam = zeros(Nblk,1);
lamMatv = [];
lamMatz = [];
for i=1:Nblk
    Ni = size(bases{i},1);
    [lam(i),ndec,lamDec] = lmivar(1,[Ni 1]);
    lamMatv = blkdiag(lamMatv, kron( eye(blk(i,1)),lamDec) );
    lamMatz = blkdiag(lamMatz, kron( eye(blk(i,2)),lamDec) );
end
lamv = lmivar(3,lamMatv);
lamz = lmivar(3,lamMatz);

% Create variables for KYP Lemma constraint
if use_kyp
    P = zeros(Nblk,1);
    for i=1:Nblk
        Nbx = order(bases{i});
        [P(i),ndec] = lmivar(1,[Nbx 1]);
    end
end

%% Create LMIs

% gamw>0
lmiterm([-1 1 1 gamwDiag],1,1);

% gamw < gampeak
lmiterm([2 1 1 gamwDiag],1,1);
lmiterm([-2 1 1 0], gampeakDiag);

% lam >0
if use_kyp
    tol = 1e-2;
    for i=1:Nblk
        [Ai,Bi,Ci,Di] = ssdata(bases{i});
        [Nbx,Nu] = size(Bi);
        lmiterm([(2+i) 1 1 P(i)],[Ai'; Bi'],[eye(Nbx) zeros(Nbx,Nu)],'s');
%         lmiterm([-(2+i) 1 1 0],[Ci'; Di']*lam(i)*[Ci Di]);
        lmiterm([-(2+i) 1 1 lam(i)], [Ci'; Di'], [Ci Di]);
        lmiterm([(2+i) 1 1 0], blkdiag(zeros(Nbx), tol * eye(Nu)));
    end
else
    for i=1:Nblk
        lmiterm([-(2+i) 1 1 lam(i)],1,1);
    end
end

% Gain Upper Bound LMI at each frequency
%  [ gI     G21    G22]  -  [ 0  0          0       ]
%  [ G21'   X        0]     [ 0  G11'X G11 G11'X G12]   > 0
%  [ G22'   0       gI]     [ 0  G12'X G11 G12'X G12]
%
% Notes: 
% 1) Let M=Mr + j Mi be a complex Hermitian matrix. 
%   a) If M=M' then Mr=Mr' and Mi = -Mi'
%   b) M>0 if and only if [Mr -Mi; Mi Mr] > 0
% 2) Let v=vr +j vi be a complex (col) vector and M=M' be a real matrix.
%   Then:
%      v' M v = (vr' - j vi') M (vr + j vi)
%             = [vr'M vr + vi' M vi] +j [vr' M vi - vi' M vr]
%             = vr'M vr + vi' M vi
%      The imaginary terms cancel because vr' M vi = vi' M vr when M=M'.
%      This relies on the fact that vr' M vi is a real scalar and hence
%      it is equal to its own transpose.
% 3) Let N be complex matrix (with more than 1 col) and M=M' be a real 
%   matrix of appropriate dimensions. Then:
%     N' M N = (Nr' - j Ni') M (Nr + j Ni)
%             = [Nr' M Nr + Ni' M Ni] +j [Nr' M Ni - Ni' M Nr]
%

G2rd = G.ResponseData(Nz+1:end,:,:);
PsivIZ = Psiv*[eye(Nv) zeros(Nv,Nd)];
PsivIZrd = PsivIZ.ResponseData;
PsizG1 = Psiz*G(1:Nz,:);
PsizG1rd = PsizG1.ResponseData;
for i=1:Nw
    % Grab data for i^th freq
    lmiid = newlmi;
    G2r = real(G2rd(:,:,i));
    G2i = imag(G2rd(:,:,i));
    PsivIZr = real(PsivIZrd(:,:,i));
    PsivIZi = imag(PsivIZrd(:,:,i));
    PsizG1r = real(PsizG1rd(:,:,i));
    PsizG1i = imag(PsizG1rd(:,:,i));
    
    % Blocks (1:2,1:2) are real blocks    
    lmiterm([-lmiid 1 1 gamwI(i)],1,1);   
    lmiterm([-lmiid 1 2 0],G2r);
    lmiterm([-lmiid 2 2 gamwZI(i)],1,1);    
    lmiterm([-lmiid 2 2 lamv],PsivIZr',PsivIZr);
    lmiterm([-lmiid 2 2 lamv],PsivIZi',PsivIZi);    
    lmiterm([lmiid 2 2 lamz],PsizG1r',PsizG1r);
    lmiterm([lmiid 2 2 lamz],PsizG1i',PsizG1i);
    
    % Blocks (3:4,3:4) are real blocks
    lmiterm([-lmiid 3 3 gamwI(i)],1,1);   
    lmiterm([-lmiid 3 4 0],G2r);
    lmiterm([-lmiid 4 4 gamwZI(i)],1,1);
    lmiterm([-lmiid 4 4 lamv],PsivIZr',PsivIZr);
    lmiterm([-lmiid 4 4 lamv],PsivIZi',PsivIZi);
    lmiterm([lmiid 4 4 lamz],PsizG1r',PsizG1r);
    lmiterm([lmiid 4 4 lamz],PsizG1i',PsizG1i);
    
    % Blocks (1:2,3:4) are imag blocks
    lmiterm([lmiid 1 4 0],G2i);
    lmiterm([-lmiid 2 3 0],G2i');
    lmiterm([-lmiid 2 4 lamz],PsizG1r',PsizG1i);
    lmiterm([lmiid 2 4 lamz],PsizG1i',PsizG1r);    
end


% Grd = G.ResponseData;
% for i=1:Nw
%     lmiid = newlmi;
%     lmiterm([-lmiid 1 1 gamwI(i)],1,1);   
%     
%     G2 = Grd(Nz+1:end,:,i);
%     lmiterm([-lmiid 1 2 0],G2);
%     
%     lmiterm([-lmiid 2 2 gamwZI(i)],1,1);
%     tmp = Psiv*[eye(Nv) zeros(Nv,Nd)];
%     lmiterm([-lmiid 2 2 lamv],tmp',tmp);
%         
%     G1 = Grd(1:Nz,:,i);
%     tmp = Psiz*G1;
%     lmiterm([lmiid 2 2 lamz],tmp',tmp);
% end

%% Create Cost Function
%  min sum_i gamw(i)
c = zeros(ndec,1);
c(1 : Nw) = 1;

%% Solve SDP
lmisys = getlmis;
opt = [0 0 0 0 1];
[copt,xopt] = mincx(lmisys,c,opt);
if isempty(copt)
	error('D-scales could not be found.');
end
gamw = xopt(1:Nw);
gamma = max(gamw);

%% Construct D-scales
Dz = [];
Dv = [];
for i=1:Nblk
    lamopt = dec2mat(lmisys,xopt,lam(i));
    Xi = bases{i}'*lamopt*bases{i};
    
    % XXX Perform spectral factorization on Xi to get Di
    [DiU,Si] = spectralfact(Xi); 
    Di = sqrt(Si)*DiU;
    Dz = blkdiag(Dz, Di * eye(blk(i, 2)));
    Dv = blkdiag(Dv, Di * eye(blk(i, 1)));
end

