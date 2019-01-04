function freq = pick_freq_grid(cl_ss, blk_str, opt, nb)
    if nargin < 3 || isempty(opt)
        freq = [];
    else
        freq = opt.FrequencyVector;
    end
    if ~isempty(freq)
        return;
    end
    if nargin < 4
        n_points = 60;
    else        
        n_points = size(blk_str(1:end - 1, :), 1) * ceil(nb * (nb + 1) / 2);
%         n_points = 3 * size(blk_str(1:end - 1, :), 1) * ceil(nb * (nb + 1) / 2);				
% 				n_points = 100;
        if n_points == 1
            n_points = 2;
        end
    end        
    Delta = reconstruct_unc(blk_str(1:end - 1, :));
    try
        wcgopt = wcOptions('VaryFrequency','on');
        [gain, ~, info] = wcgain(lft(Delta, cl_ss), wcgopt);
    catch er
%         keyboard;
%         if strcmp(er.identifier, 'MATLAB:svd:matrixWithNaNInf')
            warning(er.message);
%             save(['wcgain_grid_error_', char(datetime)], 'cl_ss', 'Delta', 'wcgopt', 'er');
            robopt = robOptions('VaryFrequency', 'on');
            [gain, ~, info] = robstab(lft(Delta, cl_ss), robopt);
%         else
%             throw(er);
%         end
    end       
% %     mubounds = mussv(cl_ss, blk_str, 's');
% %     freq_mu = mubounds.Frequency;
    freq_mu = info.Frequency;%     
    freq = freq_mu;
    freq(freq == 0) = [];
    freq(~isfinite(freq)) = [];
%     n_peaks = min(ceil(n_points / 4), numel(freq_mu));
    n_peaks = min(1, numel(freq_mu));
%     dyn_range = [freq(1) freq(end)];
% %     fc = sqrt(prod(dyn_range));
    mean_pow = mean(log10(freq));
    w_c = 10^mean_pow;
    fac = 1e4;
    w_min = max(w_c / fac, min(freq));
    w_max = min(w_c * fac, max(freq));
% %     fac = 10^4;
% %     fmin = fc / fac;
% %     fmax = fc * fac;
% %     freq(freq < fmin) = [];
% %     freq(freq > fmax) = [];
    freq = logspace(log10(w_min), log10(w_max), n_points - n_peaks)';
%     freq = sort([gain.CriticalFrequency, freq]);
    [~, ind] = sort(info.Bounds(:, 2));
    ind = flip(ind);
    freq_peak = info.Frequency(ind(1:n_peaks));
    freq_peak(~isfinite(freq_peak)) = [];
    if numel(freq_peak) ~= n_peaks
        freq_peak = [freq_peak; info.Frequency(ind(n_peaks + 1))];
    end
    freq = sort([freq; freq_peak]);
end

function Delta = reconstruct_unc(blk_str)
    Delta = [];
    for blk = 1:size(blk_str, 1)
        delta = ultidyn(['delta_dyn_', num2str(blk)], blk_str(blk, :));
        Delta = blkdiag(Delta, delta);
    end
end