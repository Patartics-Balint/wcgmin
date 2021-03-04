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
	% CommonDScale         Use a single D-scale for all the samples. This D-
	%                      scale is computed in a nonlinear optimisation
	%                      with the controller (default = 'off').
	%
	% Display              Display level (default = 'final'). Set Display to
	%                      'final' to print a one-line summary after each
	%                      outer iteration. Set Display to 'iter' to print
	%                      the iteration on the samples too. Set Display to 
	%                      'fitting' to see the printout of the D-scale
	%                      fitting as well.
	%
	% FrequencyVector      Frequency vector used for the D-scale construction
	%                      (default = []). When empty, the frequency range
	%                      and number of points are chosen automatically.
	%                      This option will set the NFrequencyPoints option
	%                      to the appropriate value.
	% 
	% HinfnormTol          Tolerance to pass when computing H-infty norm
	%                      (default = 1e-4).
	%
	% hinfstructOptions    Options to be passed to systune (default =
	%                      hinfstructOptions('RandomStart', 10,
	%                      'UseParallel', false, 'Display', 'off')).
	%
	% IterTerminateTol     Tolerance used to determinate when to terminate
	%                      the iteration (default = 0.01).
	%
	% MaxIter              Maximum number of outer iterations to perform
	%                      (default = 15).
	%
	% MaxIterDyn           Maximum number of iterations to perform on the
	%                      samples (default = 30).
	%
	% NBasisFuns           Number of basis functions used in the D-scale
	%                      construction (default = []). When empty, the total
	%                      order of the D-scales are chosen, so that it is
	%                      close to the order of the plant + controller.	
	%
	% NFrequencyPoints     Number of frequency points to be used in the
	%                      D-scale construction (default = []). The freqency
	%                      values will be chosen automatically. The
	%                      specificatin of the FrequencyVector option
	%                      overrides this setting.
	%
	% rngSeed              Sets the seed for the rng command executed at the
	%                      start of the iteration.
	%
	% systuneOptions       Options to be passed to systune (default =
	%                      systuneOptions('RandomStart', 10, 'UseParallel',
	%                      false, 'Display', 'off')).
	%
	% UseParllel           Parallel processing flag. Set to true to enable
	%                      parallel computing. This option will set the
	%                      UseParallel option of hinfstructOptions and
	%                      systuneOptions to match (default = false).
	%		
	% See also: wcgmin, hinfstructOptions, systuneOptions, hinfnorm
properties
	MaxIter = 15;
	MaxIterDyn = 30;
	Display = 'final'
	CommonDScale = 'off';
	FrequencyVector = []
	NFrequencyPoints = []
	NBasisFuns = []	
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
	
	function obj = set.FrequencyVector(obj, freq_vec)
		obj.FrequencyVector = freq_vec;
		obj.NFrequencyPoints = numel(obj.FrequencyVector);
		obj.warn_when_commonD();
	end
	
	function obj = set.NFrequencyPoints(obj, n_freq)
		if isempty(obj.FrequencyVector) || n_freq == numel(obj.FrequencyVector)
			obj.NFrequencyPoints = n_freq;
		elseif n_freq ~= numel(obj.FrequencyVector)
			obj.NFrequencyPoints = numel(obj.FrequencyVector);
			warning('FrequencyVector specified. Setting NFrequencyPoints to %d.', obj.NFrequencyPoints);
		end
		obj.warn_when_commonD();
	end
	
	function obj = set.NBasisFuns(obj, n_bf)
		if isnumeric(n_bf) && n_bf > 0 && floor(n_bf) == n_bf
			obj.NBasisFuns = n_bf;
		else
			error('The number of basis functions must be a positive integer.');
		end
		obj.warn_when_commonD();
	end
end

methods (Access = private)
	function warn_when_commonD(obj)
		if strcmp(obj.CommonDScale, 'on')
			warning('The CommonDScale option is ''on'', this option won''t take effect.');
		end
	end
end

end
