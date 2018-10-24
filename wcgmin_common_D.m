function [K_final, g_final, info] = wcgmin_Apkarian(sys, Kt, opt)
    % WCGMIN Structured robust controller design for worst-case gain optimisation.
    %
    % [K, wcg] = WCGMIN(sys, Kt) synthesises a controller K for the
    % open-loop uncertain plant sys which minimises the worst-case gain of
    % the closed loop. Kt is the tunable controller, sys is the uncertain
    % state-space model, the closed loop is lft(sys, Kt). K is the
    % controller obtained from the optimisation and wcg is the upper bound
    % on the closed-loop performance. These variables are related by:
    %       clu = lft(sys, K);
    %       g = wcgain(clu);
    %       wcg = g.UpperBound;
    %
    % [K, wcg] = WCGMIN(sys, Kt, opt) specifies options for the
    % optimisation. Use wcgminOptions to create opt.
    %
    % [K, wcg, info] = WCGMIN(sys, Kt, ...) returns an info structure
    % containing information about the algorithm execution. For the i-th
    % iteration, the info structure contains the following fields:
    %               K: Controller at the i-th iteration, ss object.
    %            gain: Worst-case gain upper bound for the closed-loop system.
    %         unc_set: Uncertainty set which is used to get samples at the i-th
    %                  iteration.
    %     info_sample: Info structure regarding the iterations on the
    %                  samples.
    %
    % See also wcgminOptions, wcgain, robstab, systune, hinfstruct, dksyn
    if nargin < 3
        opt = wcgminOptions;
    end
    rng(opt.rngSeed);
    [unc_names, dyn_unc_names, unc_nom] = get_unc(sys);
    cnt = 0;
    g_anal_ub_prev = inf;
    g_syn = nan;
    info = [];
    info_samples = [];
    unc_set = opt.UncSet;
    Drt = [];
    Dct = [];
    [~, blk_str] = components(sys, Kt);
%     cblk = blk_str(1:numel(dyn_unc_names), :);
    cblk = [];
    for name = dyn_unc_names
        name = name{1};
        cblk(end + 1, :) = size(sys.Uncertainty.(name));
    end
    nxd = ceil(order(sys) / sum(cblk(:))) + 2;
    for i = 1:size(cblk, 1)
        blk_ind = sum(cblk(1:i-1, 2)) + 1;
        d_t = tunableSS(sprintf('d_%d', i), nxd, 1, 1, 'companion');
        d_t.b.Value = zeros(nxd, 1);
        d_t.c.Value = zeros(1, nxd);
        d_t.d.Value = 1;
        Drt = blkdiag(Drt, d_t * eye(blk_str(i, 2)));
        Dct = blkdiag(Dct, d_t * eye(blk_str(i, 1)));
    end
    while true
        if cnt > 1
            g_anal_ub_prev = g_anal_ub;
        end        
        [Kt, g_anal_ub, g_anal_lb, wcunc] = analyse(sys, Kt, unc_set, unc_names, dyn_unc_names,...
            unc_nom, cnt, opt);
        print_msg(cnt, g_syn, g_anal_lb, g_anal_ub, opt);
        info = save_wcgmin_info(getValue(Kt), info_samples, g_anal_ub, unc_set, info);
        if eval_term_cond(g_syn, g_anal_lb, g_anal_ub, g_anal_ub_prev, cnt, wcunc, opt)
            g_final = min([info.gain]);
            min_ind = max(find([info.gain] == g_final));
            K_final = info(min_ind).K;
            if min_ind ~= length(info) && ~strcmp(opt.Display, 'off')
                fprintf('final: %d. wcg. anal: %.4f\n', min_ind - 1, g_final);
                info(end + 1) = info(min_ind);
            end            
            break;
        end    

        unc_set = [unc_set; wcunc];
        plant_arr = usubs(sys, unc_set);
    %     nugaps = calc_nugaps(plant_arr)
    %     sum_nugaps = sum(nugaps)
%         [K, g_syn, info_samples] = wcgmin_dyn(plant_arr, Kt, opt);

    if isa(plant_arr, 'ss')        
        [Kt, g_syn] = hinfstruct(plant_arr, Kt, opt.hinfstructOptions);        
        info_samples = struct;
    else
        [Kt, Drt, Dct, g_syn, info_samples] = synthesise(plant_arr, Kt, Drt, Dct, g_anal_ub, opt);
    end
%         Kt = setValue(Kt, K);    
        cnt = cnt + 1;    
    end
end

function [Kt, Drt, Dct, g_syn, info_samples] = synthesise(plant_arr, Kt, Drt, Dct, g, opt)
    if ~isfinite(g)
        g = 1;
    end
    g_sqrt = tunableGain('sqrt_gamma', sqrt(g));
    g_sqrt.Gain.Minimum = 0;
    g_sqrt.u = 'sqrt_gamma_u';
    g_sqrt.y = 'sqrt_gamma_y';   
    for kk = 1:nmodels(plant_arr)        
        M(:, :, kk) = components(plant_arr(:, :, kk), Kt);
    end
    clt = close_wcg_min_loop(M, Kt, Drt, Dct, g_sqrt);
    clt.u = 'clt_in';
    clt.y = 'clt_out';        
    stab_obj = TuningGoal.Gain(clt.u, clt.y, 1);
    perf_obj = TuningGoal.Gain(g_sqrt.u, g_sqrt.y, getValue(g_sqrt));
    cl = systune(append(g_sqrt, clt), perf_obj, stab_obj, opt.systuneOptions);
    Kt = setBlockValue(Kt, cl);
    Drt = setBlockValue(Drt, cl);
    Dct = setBlockValue(Dct, cl);
    [~, ~, ~, g_syn] = ssdata(getValue(g_sqrt, cl)^2);    
    info_samples = [];
end

function [Kt, g_anal_ub, g_anal_lb, wcunc] = analyse(sys, Kt, unc_set, unc_names, dyn_unc_names, unc_nom, cnt, opt)
    g_inf = 1e5;
    if cnt > 0 || opt.UseInitK
        K = getValue(Kt);
        cl_unc = lft(sys, K);
% 				[wcg, wcunc] = wcgain(cl_unc);
        [wcg, wcunc] = wcgain(cl_unc);        
        g_anal_lb = wcg.LowerBound;
        g_anal_ub = wcg.UpperBound;
        if g_anal_ub > g_inf
            [rsm, rsunc] = robstab(cl_unc);
            if rsm.UpperBound <= 1
                wcunc = rsunc;
            end            
            g_anal_ub = inf;
        end
        for j = 1:length(unc_names)
            name = unc_names{j};
            if ~isfield(wcunc, name)
                wcunc.(name) = 0;
            end
        end
        wcunc = rmfield(wcunc, dyn_unc_names);
    end
    if isempty(unc_set)
        wcunc = unc_nom;
    end    
end

function ans = eval_term_cond(g_syn, g_anal_lb, g_anal_ub, g_anal_ub_prev, cnt, wcunc, opt)
    ans = cnt >= opt.MaxIter...
        || (cnt > 0 && isempty(fieldnames(wcunc)))...
        || (is_in_interval(g_syn, [g_anal_lb, g_anal_ub], opt.IterTerminateTol)...
            && g_anal_ub >= (1 - opt.IterTerminateTol) * g_anal_ub_prev...
            && isfinite(g_anal_ub));
end

function [unc_names, dyn_unc_names, unc_nom] = get_unc(sys)
    unc = sys.Uncertainty;
    dyn_unc_names = {};
    unc_names = fieldnames(unc);
    unc_nom = struct;
    for j = 1:length(unc_names)
        name = unc_names{j};
        if isa(unc.(name), 'ultidyn')
            dyn_unc_names{end + 1} = name;
        else
            unc_nom.(name) = unc.(name).NominalValue;
        end
    end
end

function print_msg(cnt, g_syn, g_anal_lb, g_anal_ub, opt)
    if ~strcmp(opt.Display, 'off')
        if cnt == 0
            msg = 'inital ';
        else
             msg = sprintf('%d. wcg. syn.: %4.4f\t ', cnt, g_syn);
             if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
                 msg = ['\n', msg];
             end
        end    
        msg = [msg, sprintf('wcg. anal.: [%4.4f, %4.4f]\n', g_anal_lb, g_anal_ub)];
        if strcmp(opt.Display, 'iter') || strcmp(opt.Display, 'fitting')
            msg = [msg, repmat('~', 1, ceil(1.01 * length(msg))), '\n'];
        end
        fprintf(msg);
    end
end

function nugaps = calc_nugaps(plant_arr)
    n = nmodels(plant_arr);
    for k = 1:n
        for j = 1:n
            if k <= j
                Mj = plant_arr(:, :, j);
                if ~isa(Mj, 'ss')
                    Mj = lftdata(Mj);
                end
                Mk = plant_arr(:, :, k);                
                if ~isa(Mk, 'ss')
                   Mk = lftdata(Mk); 
                end                
                [~, nugap] = gapmetric(Mk, Mj);
                nugaps(k, j) = nugap;        
            end
        end
    end
    nugaps = nugaps + triu(nugaps, 1)';
end

function answer = is_in_interval(val, int, tol)
    if nargin < 3
        tol = 1e-4;
    end
    if int(1) > int(2)
        int = flip(int);
    end
    int(1) = (1 - tol) * int(1);
    int(2) = (1 + tol) * int(2);
    if val >= int(1) && val <= int(2)
        answer = true;
    else
        answer = false;
    end
end

function info = save_wcgmin_info(K, info_sample, gain, unc_set, info)
    new_info.K = K;
    new_info.info_sample = info_sample;
    new_info.gain = gain;   
    new_info.unc_set = unc_set;
    if nargin == 4
        info = new_info;
    else
        info = [info, new_info];
    end
end