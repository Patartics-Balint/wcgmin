function [K, m, info] = stabmax(sys, Kt, opt)
    if ~strcmp(opt.Display, 'off')
        fprintf('extending stability margin...\n');
    end
    if ~isa(sys, 'uss')
        error('Input system must be uss.');
    end
    rng(opt.rngSeed);
    [M, blk_str, ~, ~, nz, nw] = components(sys, Kt);
    blk_str = blk_str(1:end - 1, :);
    m_min = opt.StabilityMarginInterval(1);
    m_max = opt.StabilityMarginInterval(2);

    Dr_ss(:, :, 1:nmodels(sys)) = ss(eye(nz));
    Dc_ss(:, :, 1:nmodels(sys)) = ss(eye(nw));
    m_prev = 0;
    m_grid = 0;
    iter_cnt = 0;
    info = [];
    while true
        iter_cnt = iter_cnt + 1;
        
        K = getValue(Kt);
        cl_ss = close_stab_max_loop(M, K, Dr_ss, Dc_ss);
        m_prev = m_grid;
        for k = 1:nmodels(cl_ss)        
            if isstable(cl_ss(:, :, k))
                if hinfnorm(cl_ss(:, :, k)) ~= 0                
                    [Dr_frd{k}, Dc_frd{k}, m_grid(k), max_blk_ind] = stabmax_anal(cl_ss(:, :, k), blk_str, opt);
                else
                    m_grid(k) = inf;
                    Dr_frd{k} = eye(nz);
                    Dc_frd{k} = eye(nw);
                end
                if isfinite(m_grid(k))   
                    [Dr_ss(:, :, k), Dc_ss(:, :, k)] = fit_D_mu(Dr_frd{k}, blk_str, max_blk_ind, cl_ss(:, :, k), m_grid(k), opt);
                    m_fit(k) = 1 / hinfnorm(close_stab_max_loop(M(:, :, k), K, Dr_ss(:, :, k), Dc_ss(:, :, k)), opt.HinfnormTol);
                else
                    m_fit(k) = inf;
                    Dr_ss(:, :, k) = eye(nz);
                    Dc_ss(:, :, k) = eye(nw);                
                end
            else
                m_grid(k) = 0;
                m_fit(k) = 0;
                Dr_frd{k} = eye(nz);
                Dc_frd{k} = eye(nw);
                Dr_ss(:, :, k) = eye(nz);
                Dc_ss(:, :, k) = eye(nw);
            end
        end
        m_grid = min(m_grid(:));
        m_fit = min(m_fit(:));
        
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            fprintf('\t%d. ', iter_cnt);
            fprintf('stab. marg. anal.: %.4f (%.4f)\t', m_fit, m_grid);
        end
        
        info = save_info(K, Dr_frd, Dc_frd, Dr_ss, Dc_ss, 1/m_grid, [], info);
        
        if (m_grid >= m_min && m_fit >= m_min) || m_grid < (1 + opt.IterTerminateTol) * m_prev || iter_cnt >= opt.MaxIter
            m = m_grid;
            fprintf('\n');
            break;
        end
        
        [K, m_syn] = stabmax_syn(M, Kt, Dr_ss, Dc_ss, opt);
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            fprintf('stab. marg. syn.: %.4f\n', m_syn);
        end
        Kt = setValue(Kt, K);
    end
    if m > m_max
        [K, m, g] = reduce_stab_marg(M, Kt, Dr_ss, Dc_ss, m_max, opt);
        fprintf('\treducing stab. marg. ...\n final stab. marg.: %.4f with wcg. %.4f\n', m, g);
    end
%     print_sol_check(sys, K, m);
end

function print_sol_check(sys, K, m)
    if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting') || strcmp(opt.Display, 'final')
        stab_marg = robstab(lft(sys, K));
        fprintf('robstab: [%.4f, %.4f], design: %.4f\t', stab_marg.LowerBound, stab_marg.UpperBound, m);
        if m >= (1 - 1e-2) * stab_marg.LowerBound && m <= (1 + 1e-2) * stab_marg.UpperBound
            fprintf(':)');
        elseif m < stab_marg.LowerBound
            fprintf(':|');
        else
            fprintf(':(');
        end
        fprintf('\n');
    end
end