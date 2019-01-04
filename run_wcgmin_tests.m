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
		opt.UseParallel = true;
    try        
        [~, g, info] = wcgmin(sys, Kt, opt);
        new_results.g = g;
				new_results.info = info;
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
