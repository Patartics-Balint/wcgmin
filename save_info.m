function info = save_info(K, Dr_frd, Dc_frd, Dr_ss, Dc_ss, gain, sens, info)
    new_info.K = K;
    new_info.Dr_frd = Dr_frd;
    new_info.Dc_frd = Dc_frd;
    new_info.Dr_ss = Dr_ss;
    new_info.Dc_ss = Dc_ss;
    new_info.gain = gain;
    new_info.sens = sens;
    if nargin == 6
        info = new_info;
    else        
        info = [info, new_info];
    end    
end