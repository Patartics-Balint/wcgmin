function [Dr, Dc, gain, muub] = compute_D_skewed_mu(cl, blk, opt)
    ublk = blk(1 : end - 1, :);
    nz = sum(ublk(:, 2));
    nw = sum(ublk(:, 1));
    cl = minreal(cl, [], false);
		if nargin > 2 && ~isempty(opt.NBasisFuns)
			n_basis = opt.NBasisFuns;
		else
			n_basis = ceil(order(cl) / (nz + nw));
		end
    freqs = pick_freq_grid(cl, blk, [], n_basis);
    pp = logspace(log10(min(freqs(freqs ~= 0))), log10(freqs(end)), n_basis - 1 + 4);
    p = pp(3:end - 2);
    b = ss(1);
    for pole = p
        b = [b; zpk([], -pole, pole)];
    end
    nublk = size(blk, 1) - 1;
    bases = repmat({b}, nublk, 1);
    cl_frd = frd(cl, freqs);
    [gamma, muub, Dr, Dc] = wcgainub_lmi(cl_frd, ublk, bases, 1);
    g_check = check_smu(cl, Dr, Dc, gamma);
% 		gain = max([gamma, g_check]);
		gain = g_check;
end
