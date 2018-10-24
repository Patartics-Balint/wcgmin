function g = bisect_perf(sys, ne, nd, tol, verbose)
    % BISECT_PERF Performance gain computation using bisection.
    % 
    % Syntax
    %   g = BISECT_PERF(sys, ne, nd)
    %   g = BISECT_PERF(sys, ne, nd, tol)
    %   g = BISECT_PERF(sys, ne, nd, tol, verbose)
    %
    % The performance gain g is computed so that it is the H infinity norm of
    % sys from the last nd input channels to the last ne output channels.
    % In other words, a bisection is used to enforce
    % 
    %   ||                     +---------+                     ||
    %   || <-------------------+         +<------------------- ||
    %   ||                     |         |                     ||
    %   ||                     |   sys   |                     || = 1.
    %   ||      +---------+    |         |     +---------+     ||
    %   || <----+1/sqrt(g)<----+ ne   nd +<----+1/sqrt(g)<---- ||
    %   ||      +---------+    +---------+     +---------+     ||
    %
    if isa(sys, 'ss') && ~isstable(sys)
        g = inf;
        return;
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-4;
    end
    if nargin < 5 || isempty(verbose)
        verbose = false;
    end
    [ny, nu] = size(sys);
    nz = ny - ne;
    nw = nu - nd;
    norm_to_check = @(g)( hinfnorm(blkdiag(eye(nz), eye(ne) / sqrt(g)) * sys *...
        blkdiag(eye(nw), eye(nd) / sqrt(g)), tol) );
    g_min = hinfnorm(sys(nz+1:end, nw+1:end), tol);
    if isinf(g_min)
        g = inf;
        return;
    end
    if g_min ~= 0
        g_max = 10 * g_min;
    else
        g_max = 1;
    end
    if hinfnorm(sys(1:nz, 1:nw), tol) >= 1
        g = inf;
        return;
    end
    while norm_to_check(g_max) >= 1
        g_max = 2 * g_max;
        if g_max > 1e8 * g_min
            g = inf;
            return;
        end
    end
    while true
        g = (g_min + g_max) / 2;
        if verbose
            fprintf('\t[%d, %d]\n', g_min, g_max);
        end        
        if (g_max - g_min) <= tol * g_max  %|| isinf(g_min) || isinf(g_max)
            break;
        end
        if norm_to_check(g) >= 1
            g_min = g;
        else
            g_max = g;
        end
    end   
    g = g_max;
end

