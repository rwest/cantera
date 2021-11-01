function [class_id, job_id] = SearchID(func_name)
%% Searches for Magic # associated with a certain function name
T = readtable('FunctionDictionary/func_id.csv');
i = ismember(T{:, 1}, func_name);
ID = T{i, 2};
if isempty(ID)
    error('Method not found.');
    return;
end
class_id = floor(ID / 100);
job_id = mod(ID, 100);
end