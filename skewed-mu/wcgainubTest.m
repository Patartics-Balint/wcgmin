%% wcgainubTest
clear;
close all;

%% Create problem data
% ublk = [1 2; 3 2; 2 2];
% Ne = 2;
% Nd = 1;
% Nx = 10;
% Nb = 5;
% Nw = 10;
% blk = [ublk; Nd Ne];
% 
% rng(3);
% G = rss(Nx, Nz + Ne, Nv + Nd);
% G = stabsep(G);
% G = G * blkdiag(1 / hinfnorm(G(1:Nz, :)) * eye(Nv), eye(Nd));

load wcgmin/fitting_data_4.mat;
K = info(1).K;
usys = sys;
[M, ublk, Ne, Nd] = components(usys, K);
G = lft(M, K);
% sys = prescale(sys);
% sys = balreal(sys);
Nx = order(G);
ublk = ublk(1:end - 1, :);
blk = [ublk; Nd Ne];
% freqs = pick_freq_grid(sys, blk_str);
% freqs = [0; freqs];
Nb = 5;
w = pick_freq_grid(G, ublk);
Nw = numel(w);
% nw = 10;

Nublk = size(ublk,1);
Nv = sum( ublk(:,1) );
Nz = sum( ublk(:,2) );
if Nb==1
    bases = { ss(1) };
else
%     p = logspace(-1,1,Nb-1)';
    pp = logspace(log10(min(w(w ~= 0))), log10(w(end)), Nb - 1 + 2)';
    p = pp(2:end - 1);
    A = -diag(p);
    B = sqrt(p);
    C = diag(sqrt(p));
    bases = { [1;ss(A,B,C,0)] };
end
bases = repmat(bases, Nublk, 1);

%% Frequency grid
% w = logspace(-2, 2, Nw);
Gw = frd(G, w);

%% Solve with wcgain
[bnds, info] = mussv(Gw, blk, 's', [], [], [], 1:Nublk);
[~, ~, vlmi] = mussvextract(info);
Dz_frd = sqrtm(vlmi.Dr);
Dz_frd = Dz_frd(1:Nz, 1:Nz) / Dz_frd(end, end);
Dv_frd = sqrtm(vlmi.Dc);
Dv_frd = Dv_frd(1:Nv, 1:Nv) / Dv_frd(end, end);
muub = bnds(1,1);
muub = muub.ResponseData(:);

%% Solve with custom wcgain UB code
[gamopt, gamw, Dz, Dv] = wcgainub_lmilab(Gw,ublk,bases);

%% Display Results
figure;
subplot(3, 1, 1);
semilogx(w, muub, 'b', w, gamw, 'r--');
title('Skewed-mu');
xlabel('Frequency (rad/s)');
ylabel('UB');
legend('mussv', 'LMI', 'location', 'best');
grid on;
subplot(3, 1, 2);
Gs = blkdiag(Dz, eye(Ne) / sqrt(gamopt)) * G / blkdiag(Dv, eye(Nd) * sqrt(gamopt));
[sig_Gs, ws] = sigma(Gs);
sig_Gs = sig_Gs(1, :);
sig_Gs_w = sigma(Gs, w);
sig_Gs_w = sig_Gs_w(1, :);
semilogx(ws, sig_Gs, w, sig_Gs_w, 'r.', ws([1, end]), [1, 1], 'k');
axis tight
title('singular values of the scaled system');
xlabel('Frequency');
ylabel('sing. val.');
grid on;
subplot(3, 1, 3);
sigma(Dz_frd, 'b.-', Dz, 'r--');
xlim(Dz_frd.Frequency([1, end]));
title('D scales');
xlabel('Freq (rad/sec)');
ylabel('sing. val. (dB)');
legend('mussv', 'LMI', 'Location', 'Best');
grid on;

% close all
% semilogx(w,muub,'b',w,gamw,'r--');
% xlabel('Freq (rad/sec)');
% ylabel('UB')
% legend('mussv','LMI','Location','Best');
% grid on;
% figure();
% sigma(Dz_frd, 'b.-', Dz_ss, 'r--');
% title('D scales');
% xlabel('Freq (rad/sec)');
% ylabel('sing. val. (dB)');
% legend('mussv', 'LMI', 'Location', 'Best');
[max(muub) gamopt]