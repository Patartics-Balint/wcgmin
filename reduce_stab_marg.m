function [K, m, g] = reduce_stab_marg(M, Kt, Dr_ss, Dc_ss, m_max, opt)    
%     nz = size(Dr_ss, 1);
%     nw = size(Dc_ss, 1);
%     cl_perf_calc = close_wcg_min_loop(M, getValue(Kt), Dr_ss, Dc_ss);
%     [ny_cl, nu_cl] = size(cl_perf_calc);
%     ne = ny_cl - nz;
%     nd = nu_cl - nw;
%     for k = 1:nmodels(cl_perf_calc)
%         g(k) = bisect_perf(cl_perf_calc(:, :, k), ne, nd);
%     end
%     g = max(g);
%     if isinf(g)
%         
%     else
    g = 1;
    g_sqrt = tunableGain('g_sqrt', sqrt(g));     % Need to initialize with valid value of g? 
    [K, g] = wcgmin_dyn_syn(M, Kt, sqrt(m_max) * Dr_ss, Dc_ss / sqrt(m_max), g_sqrt, opt.systuneOptions);    
    cl = close_stab_max_loop(M, K, Dr_ss, Dc_ss);
    m_inv = hinfnorm(cl, opt.HinfnormTol);
    m = 1 ./ m_inv;
    m = min(m(:));
%     end
end