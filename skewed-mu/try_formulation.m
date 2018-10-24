clear
% close all

% ublk = [1, 2; 3, 2; 2, 2];
% ne = 2;
% nd = 1;
% nx = 10;
% nb = 5;
% nw = 10;
% 
% rng(3);
% sys = rss(nx, ne + nz, nd + nv);
% sys = stabsep(sys);
% sys = sys * blkdiag(1 / hinfnorm(sys(1:nz, :)) * eye(nv), eye(nd));

load design_data_1.mat;
% K = info(2).info_sample(15).K;
use_kyp = 1;
K = info(2).info_sample(8).K;
usys = sys;
[~, Delta] = lftdata(usys);
[M, blk, ne, nd] = components(usys, K);
sys = lft(M, K);
nx = order(sys);
ublk = blk(1:end - 1, :);
nublk = size(ublk, 1);
nz = sum(ublk(:, 2));
nv = sum(ublk(:, 1));
nb = size(ublk, 1) * ceil(order(sys) / (nz + nv)) + 1 + 2;
% nb = 3;
w = pick_freq_grid(sys, blk, [], nb);
nw = numel(w);

sysw = frd(sys, w);
[mu, muinfo] = mussv(sysw, blk, 's', [], [], [], 1:nublk);
[~, ~, vlmi] = mussvextract(muinfo);
Dz_frd = sqrtm(vlmi.Dr);
Dz_frd = Dz_frd(1:nz, 1:nz) / Dz_frd(end, end);
Dv_frd = sqrtm(vlmi.Dc);
Dv_frd = Dv_frd(1:nv, 1:nv) / Dv_frd(end, end);
muub = squeeze(mu(1).ResponseData);
wcopt = wcOptions;
wcopt.VaryFrequency = 'on';
[wcg, ~, wcinfo] = wcgain(lft(Delta, sys), w, wcopt);
wcgub = wcinfo.Bounds(:, 2);

pp = logspace(log10(min(w(w ~= 0))), log10(w(end)), nb - 1 + 4);
p = pp(3:end - 2);
b = ss(1);
for pole = p
    b = [b; zpk([], -pole, pole)];
end
bases = repmat({b}, nublk, 1);

[gamopt, gam, Dz, Dv] = wcgainub_lmilab(sysw, ublk, bases, use_kyp);
g_check = check_smu(sys, Dz, Dv, gamopt);

% [Dz, Dv, gamopt, gam] = compute_D_skewed_mu2(sys, blk);

figure;
subplot(3, 1, 1);
semilogx(w, muub, 'b.-', w, gam, 'r.--', w, wcgub, 'c:');
title('Skewed-mu');
xlabel('Frequency (rad/s)');
ylabel('UB');
legend('mussv', 'LMI', 'wcgain', 'location', 'best');
grid on;
subplot(3, 1, 2);
syss = blkdiag(Dz, eye(ne) / sqrt(gamopt)) * sys / blkdiag(Dv, eye(nd) * sqrt(gamopt));
[sig_syss, ws] = sigma(syss);
sig_syss = sig_syss(1, :);
sig_syss_w = sigma(syss, w);
sig_syss_w = sig_syss_w(1, :);
semilogx(ws, sig_syss, w, sig_syss_w, 'r.', ws([1, end]), [1, 1], 'k');
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
ylabel('sing. val.');
legend('mussv', 'LMI', 'Location', 'Best');
grid on;

[max(muub), gamopt, g_check]