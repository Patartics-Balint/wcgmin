clear
close all
rng(2);
interp_freqs = @(f)(sort([f; sqrt(f(2:end - 1) .* f(3:end))]));

% Setup system
% blk_str = [4, 4; 3, 2; 4, 2; 5, 6];
% n_blks = size(blk_str, 1);
% nwnz = sum(blk_str, 1);
% nw = nwnz(1);
% nz = nwnz(2);
% ne = 2; nd = 3;
% nx = 30;
% sys = rss(nx, nz + ne, nw + nd);
% 
% if ~isstable(sys)
%     sys = stabsep(sys);
% end
% sg = 1.01 * hinfnorm(sys(1:nz, :));
% sys = blkdiag(eye(nz) / sg, eye(ne)) * sys;
% Delta = [];
% for kk = 1:n_blks
%     Delta_blk = ultidyn(sprintf('Delta%d', kk), blk_str(kk, :));
%     Delta = blkdiag(Delta, Delta_blk);
% end
% freqs = pick_freq_grid(sys, blk_str);
% % freqs = logspace(-4, 2, 30)';
% freqs = [0; freqs(freqs ~= 0)];
% % freqs = interp_freqs(freqs);
% n_basis = 60;

% load junk2
load wcgmin/fitting_data_1.mat;
K = info(1).K;
usys = sys;
[M, blk_str, ne, nd, nz, nw] = components(usys, K);
sys = lft(M, K);
sys = prescale(sys);
sys = balreal(sys);
nx = order(sys);
blk_str = blk_str(1:end - 1, :);
n_blks = size(blk_str, 1);
freqs = pick_freq_grid(sys, blk_str);
freqs = [0; freqs];
% freqs = interp_freqs(freqs);
% freqs = [freqs; logspace(log10(freqs(end)*1.01), log10(freqs(end)) + 2, 20)'];
% freqs = logspace(-4, 5, 40);
n_basis = 100;

%% Skewed-mu computation using mussv
[mubounds, muinfo] = mussv(frd(sys, freqs), [blk_str; nd, ne], 's', [], [], [], 1:n_blks);
smu_true = mubounds(1);
[~, ~, lmidata] = mussvextract(muinfo);
gamma_true = max(mubounds(2).ResponseData);
Dr_mussv = sqrtm(lmidata.Dr);
Dc_mussv = sqrtm(lmidata.Dc);
Dr_mussv = Dr_mussv / Dr_mussv(end, end) / sqrt(gamma_true);
Dc_mussv = Dc_mussv / Dc_mussv(end, end) / sqrt(gamma_true);
Dr_mussv = Dr_mussv(1:nz, 1:nz);
Dc_mussv = Dc_mussv(1:nw, 1:nw);

opt = wcgminOptions;
opt.Display = 'fitting';
[Dr, Dc, gamma] = compute_D_skewed_mu(sys, [blk_str; nd, ne], smu_true, Dr_mussv, Dc_mussv, opt, true);