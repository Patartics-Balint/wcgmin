classdef wcgminOptions
    % WCGMINOPTIONS Creates option set for the wcgmin command.
    %
    % opt = WCGMINOPTIONS returns the default options for the wcgmin
    % command.
    %
    % opt = WCGMINOPTIONS('Option1', value1, 'Option2', value2, ...) uses
    % name/value pairs to override the default values for 'Option1',
    % 'Option2', ...
    %                                                                 
    % Supported options include:
    % 
    % Display              Display level (default = 'final'). Set Display to
    %                      'final' to print a one-line summary after each
    %                      outer iteration. Set Display to 'iter' to print
    %                      the iteration on the samples too. Set Display to 
    %                      'fitting' to see the printout of the D-scale
    %                      fitting as well.
    %
    % FrequencyVector      Frequency vector used for analysis (default =
    %                      []). When empty, the frequency range and number
    %                      of points are chosen automatically.
    % 
		% UseParllel           Parallel processing flag. Set to true to enable
		%                      parallel computing. This option will set the
		%                      UseParallel option of hinfstructOptions and
		%                      systuneOptions to match (default = false).
		%
    % HinfnormTol          Tolerance to pass when computing H-infty norm
    %                      (default = 1e-4).
    %
    % hinfstructOptions    Options to be passed to systune (default =
    %                      hinfstructOptions('RandomStart', 10,
    %                      'UseParallel', true, 'Display', 'off')).
    %
    % IterTerminateTol     Tolerance used to determinate when to terminate
    %                      the iteration (default = 0.01).
    % 
    % MaxDScalingOrder     Maximum order to be used when fitting a scalar
    %                      element of the D-scales (default = 5).
    %
    % MaxIter              Maximum number of outer iterations to perform
    %                      (default = 15).
    %
    % MaxIterDyn           Maximum number of iterations to perform on the
    %                      samples (default = 30).
    %
    % rngSeed              Sets the seed for the rng command executed at the
    %                      start of the iteration.
    %
    % systuneOptions       Options to be passed to systune (default =
    %                      systuneOptions('RandomStart', 10, 'UseParallel',
    %                      true, 'Display', 'off')).
    %
    % See also: wcgmin, hinfstructOptions, systuneOptions, hinfnorm
properties
	MaxIter = 15;
	MaxIterDyn = 30;
	Display = 'final'
	FrequencyVector = []
	DScalingOrder = []
	MaxDScalingOrder = 5
	DScalingBackoff = 1e-2
	IterTerminateTol = 1e-2
	HinfnormTol = 1e-4;
	StabilityMarginInterval = [1, 2]
	UseParallel = false
	systuneOptions = systuneOptions('RandomStart', 10, 'Display', 'off')
	hinfstructOptions = hinfstructOptions('RandomStart', 10, 'Display', 'off')
	rngSeed = 0
	UseInitK = true
	UncSet = []
end
methods
	function obj = wcgminOptions(varargin)
		cnt = 1;
		isarg = @(str)(~isempty(find(strcmp(varargin, str))));
		while cnt < nargin
			opt_name = varargin{cnt};
			try                            
				obj.(opt_name) = varargin{cnt + 1};
			catch er
				if strcmp(er.identifier, 'MATLAB:noPublicFieldForClass')
						error('Unknown oprion "%s" for the command "wcgmin".', opt_name);
				else
						throw(er);
				end
			end
			cnt = cnt + 2;
		end        
		if isarg('DScalingOrder') && isarg('MaxDScalingOrder') && ~isempty(obj.DScalingOrder)
			warning('The option MaxDScalingOrder will be ignored, because DScalingOrder was specified.');
		end
		obj.systuneOptions.UseParallel = obj.UseParallel;
		obj.hinfstructOptions.UseParallel = obj.UseParallel;
	end

	function obj = set.UseParallel(obj, use_parallel)
		obj.UseParallel = use_parallel;
		obj.systuneOptions.UseParallel = use_parallel;
		obj.hinfstructOptions.UseParallel = use_parallel;
	end
end
end