function [Kt, g] = wcgmin_dyn_syn(M, Kt, Dr_ss, Dc_ss, g_sqrt, opt)
    % WCGMIN_SYN Design the controller for the worst case minimization.
    %
    % [Kt, g] = WCGMIN_SYN(M, Kt, Dr_ss, Dc_ss, g_sqrt, opt)    
    %   Solves the optimisation probmlem below.
    %
    %    min g_sqrt  in g_sqrt, Kt
    % 
    %    subject to: g_sqrt >= 0,
    % 
    %                ||                                                   ||
    %                ||       +----+      +-----------+    +---------+    ||
    %                || <-----+ Dr <------+           <----+ inv(Dc) <--- ||
    %                ||       +----+      |           |    +---------+    ||
    %                ||                   |           |                   ||
    %                ||    +----------+   |     M     |    +----------+   ||
    %                || <--+ 1/g_sqrt <---+           <----+ 1/g_sqrt <-- ||
    %                ||    +----------+   |           |    +----------+   || <= 1
    %                ||                   |           |                   ||
    %                ||          +--------+           <--------+          ||
    %                ||          |        +-----------+        |          ||
    %                ||          |                             |          ||
    %                ||          |            +----+           |          ||
    %                ||          +------------> Kt +-----------+          ||
    %                ||                       +----+                      ||
    %                ||                                                   ||    
    %
    % [K, g] = WCGMIN_SYN(M, ny, nu, Dr_ss, Dc_ss, opt)
    %   NOT IMPLEMENTED YET
    if nargin == 6
        g_sqrt.Gain.Minimum = 0;
        g_sqrt.y = 'g_sqrt_out';
        g_sqrt.u = 'g_sqrt_in';
        clt = close_wcg_min_loop(M, Kt, Dr_ss, Dc_ss, g_sqrt);
        clt.u = 'clt_in';
        clt.y = 'clt_out';        
        stab_obj = TuningGoal.Gain(clt.u, clt.y, 1);
        perf_obj = TuningGoal.Gain(g_sqrt.u, g_sqrt.y, getValue(g_sqrt));
        cl = systune(append(g_sqrt, clt), perf_obj, stab_obj, opt);
%         K = getValue(Kt, cl);
        Kt = setBlockValue(Kt, cl);
        [~, ~, ~, g] = ssdata(getValue(g_sqrt, cl)^2);
    end
    if nargin == 4
    end
end