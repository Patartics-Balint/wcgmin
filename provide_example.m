function [sys, ny, nu, nxK] = provide_example(type, ex_no)
	% Provide example system for sutructured control desing against mixed
	% uncertainy.
	%
	% Usage
	%   [usys, ny, nu, nxK] = PROVIDE_EXAMPLE(type, ex_no)
	% Inputs:
	%   type: 'collected' or 'random'
	%   ex_no: if type == 'collected' ex_no is the ordinal number of the
	%   example, if type == 'random' ex_no is the seed for the random number
	%   generation
	% Outputs:
	%   usys: uncertain system
	%   ny: number of measured outputs
	%   nu: number of control inputs
	%   nxK: order of the structured controller
	% See also wcgmin
	
    if nargin < 2
        if strcmp(type, 'collected')
            ex_no = 12;
        end
        if strcmp(type, 'random')
            ex_no = 0;
        end
    end
    switch type
        case 'random'
					[M, Delta, ny, nu, nxK] = provide_random_examle(ex_no);
        case 'collected'
					[M, Delta, ny, nu, nxK] = provide_collected_examle(ex_no);			
    end
    sys = lft(Delta, M);
end

function [M, Delta, ny, nu, nxK] = provide_random_examle(ex_no)
	nxK = 2;
	blk_str = [2, 2; 1, 1; 3, 3];
	Delta = [];
	for i = 1:size(blk_str, 1)
			Deltai = ultidyn(['Delta_' int2str(i)], blk_str(i, :));
			Delta = append(Delta, Deltai);
	end
	[nw, nz] = size(Delta);
	ne = 2; nd = 1;
	ny = 3; nu = 2;
	nx = 6;
	rng(ex_no);
	M = rss(nx, nz + ne + ny, nw + nd + nu);
	%     M = canon(M, 'modal');
	%     M.a = blkdiag(eye(nx - 2), -1 * eye(2)) * M.a;
	%         M = blkdiag(0.035 * eye(nz), eye(ne), eye(ny)) * M;
	M = blkdiag(0.11 * eye(nz), 0.1 * eye(ne), eye(ny)) * M;
	%         M = blkdiag(0.1 * eye(nz), 0.1 * eye(ne), eye(ny)) * M;
end

function [M, Delta, ny, nu, nxK] = provide_collected_examle(ex_no)
	nxKs = [2, 1, 4, 3, 5, 1, 4, 2, 2, 2, 1, 3, 1, 2, 3, 4, 4, 3, 1, 2, 1, 3, 4, 4, 3, 1, 3, 4, 2, 5, 1];
	nxK = nxKs(ex_no);
	load(['example', num2str(ex_no)]);
	Nblk = size(blk,1);
	Delta = [];
	Nr = 0;
	Nd = 0;
	for i=1:Nblk
			blki = blk(i,:);
			if blki(1)<0
					% UREAL
					Nr = Nr+1;
					I = eye( abs(blki(1,1)) );
					Deltai = ureal(['delta_real_' int2str(Nr)],0)*I;
			else
					% ULTIDYN
					Nd = Nd +1;
					Deltai = ultidyn(['delta_dyn_' int2str(Nd)],blki);
			end
			Delta = blkdiag(Delta,Deltai);
	end
	M = P;  
end