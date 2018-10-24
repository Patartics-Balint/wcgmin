clear;
close all;

rng(0);
ne = 2;
nd = 2;
ny = 2;
nu = 2;
nx = 7;

blk_str = [1, 1; 2, 2; 1, 1; nd, ne];
nw = sum(blk_str(1:end-1, 1));
nz = sum(blk_str(1:end-1, 2));
plant = rss(nx, nz + ne + ny, nw + nd + nu);
plant = blkdiag(0.2 * eye(nz), 0.2 * eye(ne), eye(ny)) * plant;
n_unc_blks = size(blk_str(1:end - 1, :), 1);
delta = [];
for kk = 1:n_unc_blks
    delta = blkdiag(delta, ultidyn(sprintf('delta%d', kk), blk_str(kk, :)));
end

[K, ~, gamma] = hinfsyn(plant, ny, nu);

cl = lft(plant, K);
freq = logspace(-2, 3, 100);
cl_frd = frd(cl, freq);
mubounds = mussv(cl_frd, blk_str, 's', [], [], [], 1:n_unc_blks);

mu_mussv = squeeze(mubounds(2).ResponseData)';

mu_lmi = [];
D = [];
for kk = 1:numel(freq)
    G = cl_frd.ResponseData(:, :, kk);
    g = sdpvar(1);
    X = [];
    for ii = 1:n_unc_blks
       X = blkdiag(X, sdpvar(1) * eye(blk_str(ii, 1)));
    end

    lmis_lhs = {};
    lmis_lhs{end + 1} = -X;
    lmis_lhs{end + 1} = -g;
    lmis_lhs{end + 1} = G' * blkdiag(X, g*eye(ne)) * G - blkdiag(X, eye(nd));
    lmis = [];
    for ii = 1:numel(lmis_lhs)
        lmis = [lmis, lmis_lhs{ii} <= 0];
    end
    opt = sdpsettings;
    opt.solver = 'lmilab';
    res = optimize(lmis, -g, opt);
    D(:, :, kk) = sqrtm(value(X));
    mu_lmi(kk) = 1/sqrt(value(g));
end
semilogx(freq, mu_mussv, freq, mu_lmi, '--', 'LineWidth', 2);
xlabel('frequency [rad/s]');
ylabel('skewed-mu');