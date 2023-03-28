function [Dr_frd, Dc_frd, g, sens, smu] = wcgmin_dyn_anal(cl_ss, blk_str, opt)
    % WCGMIN_ANAL Perform the analysis in the worst case minimisation.
    %
    % [Dr_frd, Dc_frd, g, sens] = WCGMIN_ANAL(cl, blk_str, freq)
    %   Solves the skewed-mu optimisation problem depicted below on the
    %   frequencies in freq.
    %
    %    min g_sqrt  in g_sqrt, Kt
    % 
    %    subject to: g_sqrt >= 0,
    % 
    %                ||                                                       ||
    %                ||        +-----+      +-----------+    +----------+     ||
    %                ||  <-----+ Dru <------+           <----+ inv(Dcu) <---  ||
    %                ||        +-----+      |           |    +----------+     ||
    %                ||                     |           |                     ||
    %                ||    +------------+   |     M     |    +------------+   ||
    %                || <--+ Drp/g_sqrt <---+           <----+ Dcp/g_sqrt <-- ||
    %                ||    +------------+   |           |    +------------+   || <= 1
    %                ||                     |           |                     ||
    %                ||            +--------+           <--------+            ||
    %                ||            |        +-----------+        |            ||
    %                ||            |                             |            ||
    %                ||            |            +----+           |            ||
    %                ||            +------------> Kt +-----------+            ||
    %                ||                         +----+                        ||
    %                ||                                                       || 
    %
    %   * Dr_frd and Dc_frd are obtainded from the scales resulting from
    %     the optimisation. Drp = drp * I, Dcp = dcp * I and Dr_frd = Dru /
    %     drp, Dc_frd = Dru \ dcp.
    %   * g is the gain resulting from the optimisation.
		cl_ss = minreal(cl_ss, [], false);
    freq = pick_freq_grid(cl_ss, blk_str, opt);        
    cl_frd = frd(cl_ss, freq);
    n_unc_blks = size(blk_str, 1) - 1;
		[mubounds, muinfo] = mussv_wrapper(cl_frd, blk_str, 's', 1:n_unc_blks);
    smu = mubounds(1);
    g = max(mubounds(1).ResponseData);
    [~, ~, lmiscales] = mussvextract(muinfo);
    nw = sum(blk_str(1:end-1, 1));
    nz = sum(blk_str(1:end-1, 2));
    Dr_frd = get_D_unc(lmiscales.Dr, g, nz);
    Dc_frd = get_D_unc(lmiscales.Dc, g, nw);
    muinfo_sens = muinfo.sens;
    s = muinfo_sens.ResponseData;
    muinfo_sens_max = max(s(:));
    muinfo_sens = muinfo_sens/muinfo_sens_max;
    sens = muinfo_sens;
%     sens1 = calc_sens(cl_ss, Dr_frd, Dc_frd, blk_str(1:end-1, :), g);    
%     sens3(:, 1:size(blk_str, 1)) = mubounds(1) / g;   
%     sens = sens3;
%     sens = sens2.*sens3;
%     sens = sens1;    
    norm_frd = hinfnorm(close_wcg_min_loop(cl_frd, [], Dr_frd, Dc_frd, sqrt(g)), opt.HinfnormTol);
    if norm_frd > 1 || norm_frd < 0.99
        error('The soulution to the mu-analysis is incorrect.');
    end
end

function D_unc = get_D_unc(D_squared, g_mu, n_unc_ch)
    D = sqrtm(D_squared);
    D = D / D(end, end) / sqrt(g_mu);
    D_unc = D(1:n_unc_ch, 1:n_unc_ch);
end