function [Dr_frd, Dc_frd, m, max_blk_ind] = stabmax_anal(cl_ss, blk_str, opt)
    % STABMAX_ANAL Design the controller for the stability maximisation.
    %
    % [Dr_frd, Dc_frd, m, max_blk_ind] = STABMAX_ANAL(cl, blk_str, opt)
    %   Solves the optimisation problem below on a frequency grid.
    % 
    %                 ||                                          ||
    %                 ||     +----+   +--------+   +---------+    ||
    %                 ||     |    |   |        |   |         |    ||
    %     1/m =  min  || <---+ Dr +---+   cl   +---+ inv(Dc) +--+ ||
    %                 ||     |    |   |        |   |         |    ||
    %          Dr, Dc ||     +----+   +--------+   +---------+    ||
    %                 ||                                          ||
    %
    %   * Dr_frd and Dc_frd are the scales on a frequency grid after the
    %     division by the blockdiagonal element corresponding to max_bld_ind
    %   * m is the inverse of the H-infinty norm of the saled system (i.e.
    %     the stability margin of cl)
    %   * max_blk_ind is the index of the uncertainty block
    %     blk_str(max_bld_ind, :) for which the scales are the biggest.
%     if isempty(opt.FrequencyVector)
%         mubounds = mussv(cl_ss, blk_str, 's');
%         freq = mubounds.Frequency;
%     else
%         freq = opt.FrequencyVector;
%     end
    freq = pick_freq_grid(cl_ss, blk_str, opt);
    [mubounds, muinfo] = mussv(frd(cl_ss, freq), blk_str, 's');
    [~, ~, lmiscales] = mussvextract(muinfo);
    [~, max_blk_ind] = max(sum(blk_str, 2));
    Dr_frd = get_D_frd(lmiscales.Dr, blk_str(:, 2), max_blk_ind);
    Dc_frd = get_D_frd(lmiscales.Dc, blk_str(:, 1), max_blk_ind);    
    m = 1 / max(mubounds(1).ResponseData);
end

function D_frd = get_D_frd(D_squared, blk_str_col, max_blk_ind)
    D_frd = sqrtm(D_squared);
    max_blk_elem_ind = sum(blk_str_col(1:max_blk_ind));
    d_max_blk = D_frd(max_blk_elem_ind, max_blk_elem_ind);
    D_frd = D_frd / d_max_blk;
end