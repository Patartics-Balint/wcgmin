function [Dr_ss, Dc_ss, nxD] = fit_D_mu(Dr_frd, blk_str, max_blk_ind, cl, m_frd, opt)
    if ~isempty(opt.DScalingOrder)
        orders = opt.DScalingOrder;
    else
        orders = 0:opt.MaxDScalingOrder;
    end    
    for nxD = orders    
        Dr_ss = [];
        Dc_ss = [];
        for blk_ind = 1:size(blk_str, 1)
            blk_size_r = blk_str(blk_ind, 2);
            blk_size_c = blk_str(blk_ind, 1);            
            if blk_ind == max_blk_ind
                d_ss = 1;
            else
                blk_start_ind = sum(blk_str(1:blk_ind - 1, 2)) + 1;
                d_frd = Dr_frd(blk_start_ind, blk_start_ind);
                d_ss = fitfrd(genphase(d_frd), nxD, 0, [], 1);
            end
            Dr_ss = blkdiag(Dr_ss, d_ss * eye(blk_size_r));
            Dc_ss = blkdiag(Dc_ss, d_ss * eye(blk_size_c));
        end
        m_fit = 1 / hinfnorm(close_stab_max_loop(cl, [], Dr_ss, Dc_ss));
        degr = (m_frd - m_fit) / m_fit;
        if strcmp(opt.Display, 'fitting')
            fprintf('\t\tD order: %d\t stab. marg. degradation: %.2f%%\n', nxD, degr * 100);
        end
        if abs(degr) <= opt.DScalingBackoff
            break;
        end
    end
end