function [K, g, m, nxK] = wcg_dksyn(sys, ny, nu, tol, verbose)
    % WCG_DKSYN Perform D-K synthesys with bisection for worst-case gain
    % minimisation
    %
    % [K, g, m, nxK] = WCG_DKSYN(sys, ny, nu, tol, verbose)
    if nargin < 4 || isempty(tol)
        tol = 1e-2;
    end
    if nargin < 5
        verbose = true;
    end
    lower_bound = 1e-2;
    upper_bound = 1e2;
    fact = 1e2;
    while true
        [K, g, rob_perf, nxK] = bisect_wcgain(sys, ny, nu, tol, verbose, lower_bound, upper_bound);
        if g - lower_bound < tol
            upper_bound = lower_bound;
            lower_bound = lower_bound / fact;            
        elseif upper_bound - g < tol
            lower_bound = upper_bound;
            upper_bound = upper_bound * fact;
        else
            break;
        end
        if lower_bound < 1e-8 || upper_bound > 1e8
            warning('The gain is beyond reasonable limits. Stopping search.');
            break;
        end        
    end
    m = 1 \ rob_perf;    
    if verbose
        cl = lft(sys, K);
        sm = robstab(cl);
        wcg = wcgain(cl);
        fprintf('robsatb: [%.4f, %.4f], design: %.4f\n', sm.LowerBound, sm.UpperBound, m);
        fprintf('wcgain: [%.4f, %.4f], design: %.4f\n', wcg.LowerBound, wcg.UpperBound, g);
    end    
end

function [K, g, rob_perf, nxK] = bisect_wcgain(sys, ny, nu, tol, verbose, lb, ub)
    [ny_sys, nu_sys] = size(sys);
    ne = ny_sys - ny;
    nd = nu_sys - nu;
    K = zeros(nu, ny);
    g_min = lb;
    g_max = ub;
    while true
        g = (g_min + g_max) / 2;
        P = blkdiag(eye(ne) / sqrt(g), eye(ny)) * sys * blkdiag(eye(nd) / sqrt(g), eye(nu));
        [K, ~, rob_perf] = dksyn(P, ny, nu);
        nxK = length(K.a);
        if verbose
            fprintf(' [%d, %d]\trob. perf.: %.4f\tnxK = %d\n', g_min, g_max, rob_perf, nxK);
        end
        if (g_max - g_min) <= tol
            break;
        end
        if rob_perf >= 1
            g_min = g;
        else
            g_max = g;
        end
    end
    g = g_max;
end