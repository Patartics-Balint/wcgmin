% Run wcgmin on the system in the examples folder. The resulting printout
% and data is saved to ./'test reports' for later comparison.

clc
clear
close all

test_path = which(mfilename);
path_parts = strsplit(test_path, '/');
test_dir = strjoin([path_parts(1 : end - 1), 'test_reports'], '/');
cd(test_dir);

warning('on');
warning('backtrace', 'off'); % file links screw up the report
date = strrep(char(datetime), ' ', '_');
report_filename = ['test_report_', date, '.txt'];
diary(report_filename);
results = [];
examples = 1:31;
tmax = 10 * 60;
datafile = ['results_', date, '.mat'];
if exist(datafile)
    load(datafile);
end
for example = examples
	if strcmp(get(0, 'Diary'), 'off')
		diary(report_filename);
	end
	if example > examples(1)
		fprintf('\n\n');
	end
	fprintf('## Test case with example no. %d ############################\n', example);    
	[sys, ny, nu, nxK] = provide_example('collected', example);
	new_results.ex_no = example;
	Kt = tunableSS('K', nxK, nu, ny, 'companion');
	opt = wcgminOptions;
	opt.Display = 'iter';
	opt.systuneOptions.Display = 'final';
	opt.hinfstructOptions.Display = 'final';
	opt.UseParallel = true;
	opt.CommonDScale = 'on';
	try
		tic;
			[~, g, info] = wcgmin(sys, Kt, opt);
% 			[~, ~, g, info] = fetchNext(parfeval(@wcgmin, 3, sys, Kt, opt), tmax);
% 			if isempty(g) && isempty(info)
% 				error('The example takes more than %d seconds to compute.', tmax);
% 			end
		t_syn =	ceil(toc);
		new_results.g = g;
		new_results.info = info;
		new_results.t_syn = t_syn;
	catch er
		fprintf(1, 'ERROR: %s\n%s\n', er.identifier, er.message);
		new_results.g = nan;               
	end
	new_results.nxK = nxK;    
	if isempty(results)
		results = new_results;
	else
		results(example) = new_results;
	end
	save(datafile, 'results');
	diary('off');
end
fprintf('\n\n');
res_table = struct2table(results);
diary('off');
