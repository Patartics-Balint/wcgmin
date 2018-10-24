function compare_results(new_datafile, old_datafile)
    new_date = extract_date(new_datafile);
    old_date = extract_date(old_datafile);
    load(new_datafile);
    new_res = results;
    load(old_datafile);
    old_res = results;
    clear('results');
    available_ex = min(length(new_res), length(old_res));
    fprintf('Ex. \t %s -> %s\n', old_date, new_date);
    for k = 1:available_ex
        [diff_num, diff_op] = eval_diff(new_res(k), old_res(k));
        fprintf('%d \t %s \t', new_res(k).ex_no, diff_op);
        if length(diff_op) < 6
            fprintf('\t');
        end
        fprintf(' (%s)\n', diff_num);
    end
end

function date = extract_date(datafile)
    components = strsplit(datafile, '_');
    date = components{2};
end

function [diff_num, diff_op] = eval_diff(new_res, old_res)
    same_tol = 0.05;
    much_tol = 1;
    g_new = new_res.g;
    g_old = old_res.g;
    if isfinite(g_new) && isfinite(g_old)
        diff = (g_new - g_old) / g_old;
        if abs(diff) <= same_tol
            diff_op = 'same';
        elseif diff > 0
            diff_op = 'worse';
        elseif diff < 0
            diff_op = 'better';
        end
        if abs(diff) > much_tol
            diff_op = ['much ', diff_op];
        end
        diff_num = sprintf('%.2f%%', diff*100);
    elseif isfinite(g_new) && ~isfinite(g_old)
        diff_op = 'much better';
        diff_num = 'bacame finite';
    elseif ~isfinite(g_new) && isfinite(g_old)
        diff_op = 'much worse';
        diff_num = 'became infinite';
    else
        diff_op = 'same';
        diff_num = 'infinite';
    end   
end