%% This file recreates the examples in Patartics, Seiler, Vanek: Structured
% Robust Synthesis with Parameter-Dependent D-Scales. Submitted to the
% American Control Conference, 2019. The results in the paper were optained
% using 64-bit Ubuntu 16.04, Intel® Core™ i5-4590 CPU @ 3.30GHz × 4 with
% 11,6 GiByte memory. There are random processes involved in the
% algorithms, therefore the results of this script may differ slightly form
% the numbers published.
close all
clear

%% Setup of the benchmark example in Section IV-A

[P_A, ny, nu, nxK] = provide_example('collected', 12);

%% Control desing for the example in Section IV-A
% The design algorithm takes a considerable amount of time to finish.

% tunable controller with 'companion' parametrisation
Kt_A = tunableSS('K', nxK, nu, ny, 'companion');
% options for the desing algorithm
opt = wcgminOptions;
% controller desing using the algorithm
[K_A, g_A, info_A] = wcgmin(P_A, Kt_A, opt);

%% Setup of the example in Secton IV-B

% nominal plant
G = tf(1, [1, 1]);
% paramteric and dynamic uncertainties
delta_r = ureal('delta_r', 0, 'Range', [-1, 1]);
delta_c = ultidyn('delta_c', [1, 1]);
% formation of the uncertain plant
Wdel = 0.05;
Del_tilde = 1/(1 + 0.98 * delta_r) * sqrt(Wdel) * delta_c * sqrt(Wdel) *...
	1/(1 - 0.98 * delta_r);
Gu = G + Del_tilde;
% weighting functions for the generalised plant
We = tf(1, [1/10, 1]);
Wd = tf(1, [1/10, 1]);
Wn = 1e-3 * makeweight(0.99, 10, 1e3);
Wu = 0.1 * makeweight(0.99, 10, 100);
% formation of the generalised plant
P = [[0, 0, Wu]; [0; Wn], [We; 1] * Gu * [Wd, 1]];

%% Analysis of the example in Section IV-B

% plot worst-case gain of Gu over the dynamic uncertainty as a function of 
% the paramteric uncertainty
del_r = linspace(-1, 1, 101); % paramteric uncertainy samples
wcg_samples = usubs(Gu, 'delta_r', del_r);
wcg = [];
for kk = 1 : nmodels(wcg_samples)
	% compute worst-case gain for a sample
	wcg_temp = wcgainWithBNB(wcg_samples(:, :, kk));	
	wcg(end + 1) = wcg_temp.UpperBound;
end
figure;
plot(del_r, wcg);
xlabel('delta_r');
ylabel('worsd-case gain');

% plot D-scales corresponding to the worst-case gain upper bound of Gu for
% a few values of the paramteric uncertainy
del_r = linspace(-1, 1, 5);
d_samples = usubs(Gu, 'delta_r', del_r);
% block structure of the dynamic uncertainty and the performance block
blk_str = [1, 1; 1, 1];
% frequency grid for the computation
freq = logspace(-1, 2, 100);
% pull out the uncertainty into an LFT
T = lftdata(d_samples);
% Call mussv for the skewed-mu computation. The last arguments are
% undocumented.
[mubounds, muinfo] = mussv_wrapper(frd(T, freq), blk_str, 's', 1);
% get the D-scales
[~, ~, vlmi] = mussvextract(muinfo);
% transform the D-scales to the form they appear in the paper
Dd = sqrtm(vlmi.Dr);  
Dd = Dd(1, 1) / Dd(end, end);
% plot magnitudes
figure;
sigma(Dd);

%% Desing controller for the example in Section IV-B
% The design algorithm takes a considerable amount of time to finish.

% tunable controller initialised to be 0
Kt_A = tunableSS('Kt', ss(-1, 0, 0, 0));
% options for the algorithm
opt = wcgminOptions;
% controller desing using the algorithm
[K_B, g_B, info_B] = wcgmin(P, Kt_A, opt);

%% Analysis of the closed loop in Section IV-B

nxD = 4; % state order of the D-scales
% tunable D-scale for the computation
Ddt = tunableSS('Dd', nxD, 1, 1);
% initialise the D-scale to identity
Ddt = setValue(Ddt, ss(Ddt.a.Value, zeros(nxD, 1), zeros(1, nxD), 1));
% tunable performance gain
g_init = 1;
g_sqrt = tunableGain('sqrt_gamma', sqrt(g_init));
g_sqrt.Gain.Minimum = 0;
g_sqrt.u = 'sqrt_gamma_u';
g_sqrt.y = 'sqrt_gamma_y';

% parametric uncertainty samples
Ps = usubs(P, 'delta_r', [0, 1]); % change [0, 1] to [-1, 1] to get the
% other result in the paper
Ms = lftdata(Ps);
% form the constratint of the optimisation
scaled_sys = blkdiag(Ddt, 1 / g_sqrt * eye(2)) * lft(Ms, K_B) *...
	blkdiag(inv(Ddt), 1 / g_sqrt * eye(2));
scaled_sys.u = 'constratint_u';
scaled_sys.y = 'constratint_y';
constratint = TuningGoal.Gain(scaled_sys.u, scaled_sys.y, 1);
perf_obj = TuningGoal.Gain(g_sqrt.u, g_sqrt.y, getValue(g_sqrt));
% run optimisation
cl = systune(append(g_sqrt, scaled_sys), perf_obj, constratint,...
	opt.systuneOptions);
% get the value of the obtained gain
[~, ~, ~, g_common_d] = ssdata(getValue(g_sqrt, cl)^2);
% get the common D-scale
Dd_common = getValue(Ddt, cl);