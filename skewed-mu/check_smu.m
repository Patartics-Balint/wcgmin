function smu = check_smu(cl, Dz, Dw, smu_ref)
    sqrt_gamma = tunableGain('sqrt_gamma', sqrt(smu_ref));
    sqrt_gamma.Gain.Minimum = 0;
    sqrt_gamma.u = 'sqrt_gamma_u';
    sqrt_gamma.y = 'sqrt_gamma_y';
    clt = close_wcg_min_loop(cl, [], Dz, Dw, sqrt_gamma);
    clt.u = 'clt u';
    clt.y = 'clt y';
    soft = TuningGoal.Gain(sqrt_gamma.u, sqrt_gamma.y, smu_ref);
    hard = TuningGoal.Gain(clt.u, clt.y, 1);
    systune_opt = systuneOptions;
    systune_opt.Display = 'off';
    systune_opt.UseParallel = false;
%     systune_opt.RandomStart = 10;
    cl = systune(append(sqrt_gamma, clt), soft, hard, systune_opt);
    [~, ~, ~, smu] = ssdata(getValue(sqrt_gamma, cl)^2);    
end