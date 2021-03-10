function [K_final, g_final, info_wcg] = wcgmin_dyn(sys, Kt, opt)
    if nargin < 3
        opt = wcgminOptions;
    end
    rng(opt.rngSeed);
    if isa(sys, 'ss')        
        [Kt, g_final] = hinfstruct(sys, Kt, opt.hinfstructOptions);
        K_final = Kt;
        info_wcg = struct;
        return;
    end
    g_sqrt = tunableGain('g_sqrt', sqrt(1));
    g_prev = inf;
    g_grid = g_prev;
    iter_cnt = 1;
    info_wcg = [];                        
    while true        
        g_prev = g_grid;
        [Dr_frd, Dc_frd, Dr_ss, Dc_ss, g_fit, g_grid, sens] = analyse(sys, Kt, iter_cnt, opt);
        info_wcg = save_info(Kt, Dr_frd, Dc_frd, Dr_ss, Dc_ss, g_grid, sens, info_wcg);        
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            fprintf('%d. ', iter_cnt);
            fprintf('wcg. anal.: %.4f (%.4f)\t', g_fit, g_grid);        
        end        
        if (g_grid >= (1 - opt.IterTerminateTol) * g_prev && isfinite(g_grid))...
                || iter_cnt > opt.MaxIterDyn
            g_final = min([info_wcg.gain]);
            min_ind = max(find([info_wcg.gain] == g_final));
            K_final = info_wcg(min_ind).K;
            if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
                fprintf('\n');
            end
            if min_ind ~= length(info_wcg) &&...
                (strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting'))
                fprintf('final: %d. wcg. anal: %.4f\n', min_ind, g_final);
            end
            break;
        end 
        [Kt, g_syn] = synthesise(sys, Kt, Dr_ss, Dc_ss, g_sqrt, g_fit, opt);
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
                fprintf('wcg. syn.: %.4f\n', g_syn);
        end        
        iter_cnt = iter_cnt + 1;       
    end    
end

function [Dr_frd, Dc_frd, Dr_ss, Dc_ss, g_fit, g_grid, sens] = analyse(sys, Kt, iter_cnt, opt)
    [M, blk_str, ~, ~, nz, nw] = components(sys, Kt);
    K = getValue(Kt);
    cl = lft(M, K);
    try
			if opt.UseParallel
				parfor k = 1:nmodels(cl)
            [Dr_frd{k}, Dc_frd{k}, g_grid(k), sens{k}, smu] =...
							wcgmin_dyn_anal(cl(:, :, k), blk_str, opt);                
            [Dr_ss(:, :, k), Dc_ss(:, :, k), g_fit(k)] =...
							compute_D_skewed_mu(cl(:, :, k), blk_str, g_grid(k), opt);
				end
			else
        for k = 1:nmodels(cl)
            [Dr_frd{k}, Dc_frd{k}, g_grid(k), sens{k}, smu] =...
							wcgmin_dyn_anal(cl(:, :, k), blk_str, opt);                
            [Dr_ss(:, :, k), Dc_ss(:, :, k), g_fit(k)] =...
							compute_D_skewed_mu(cl(:, :, k), blk_str, g_grid(k), opt);
				end
			end
    catch er_anal
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            fprintf(['analysis attempt failed: ', er_anal.message, '\n']);
        end
        if iter_cnt > 2
            rethrow(er_anal);
        end
        if ~exist('g_syn')
            Dr_ss = repmat(eye(nz), 1, 1, nmodels(M));
            Dc_ss = repmat(eye(nw), 1, 1, nmodels(M));
            Dr_frd = Dr_ss;
            Dc_frd = Dc_ss;
            g_fit = inf;
            g_grid = g_fit;
            sens = [];
        else
            [K, m, ~] = stabmax(sys, Kt, opt);
            if m < 1
                error('The stability margin is %.4f < 1, worst-case gain minimisation is impossible.', m);
            end
            Kt = setValue(Kt, K);                
        end
    end        
    g_grid = max(g_grid(:));
    g_fit = max(g_fit(:));
end

function [Kt, g_syn] = synthesise(sys, Kt, Dr_ss, Dc_ss, g_sqrt, g_fit, opt)
    M = components(sys, Kt);
    if isfinite(g_fit)
       g_init = g_fit;
    else
       g_init = 1;
    end
    g_sqrt = setValue(g_sqrt, sqrt(g_init));     
    try
       [Kt, g_syn] = wcgmin_dyn_syn(M, Kt, Dr_ss, Dc_ss, g_sqrt, opt.systuneOptions);
    catch er_syn
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            fprintf(['sysnthesys attempt failed: ', er_syn.message, '\n']);
        end
        if iter_cnt > 2
            rethrow(er_syn);
        end            
        g_syn = inf;
        [K, m] = stabmax(sys, Kt, opt);
        if m < 1
            error('The stability margin is %.4f < 1, worst-case gain minimisation is impossible.', m);
        end            
    end
%     Kt = setValue(Kt, K);    
end