%% Script that runs the petab_select tests and evaluates the results
fprintf(2, 'TEST FOR PETAB EXTENSION >>PETAB-SELECT<<\n');

% define test case folders
cases = 1:9;
nCases = length(cases);
names = arrayfun(@(x) sprintf('%04d', x), cases, 'UniformOutput', false);
mainDir = pwd();

% create results data structures
results = createEmptyResultsTable(names);


%% Run test cases and evaluate individual results
try
    parpool();
catch
end

for idx = 1:nCases
    name = names{idx};
    cd(fullfile(mainDir, name));
    try
        runExample;
        results = evaluateTestResult(results, idx);
        if results{idx, "passed"}==0
             fprintf('test case %s failed\n', name)
        end
    catch ME
        fprintf('error during test of case %s\n', name)
        results{idx, "passed"} = -1;
        results{idx, "error"} = string(getReport(ME, "extended", "hyperlinks", "off"));
    end
    cd(mainDir)
end
writetable(results, fullfile(mainDir, 'TestResults.csv'));

%% Summarize results 

if all(results.passed==1)
    fprintf(2, 'Correct results all test cases: %s\n', strjoin(names(results.passed==1)));
    fprintf(2, 'PASSED\n' );
else    
    fprintf(2, 'Correct results in test case(s) %s\n', strjoin(names(results.passed==1)));
    fprintf(2, 'Results not identical in test case(s) %s\n', strjoin(names(results.passed==0)));
    fprintf(2, 'Errors in test case(s) %s\n', strjoin(names(results.passed==-1)));
    error('FAILED');
end


%% SUBFUNCTIONS

% -------------------------------------------------------------------------
function results = createEmptyResultsTable(names)
% create an empty table with the results of the tests

emptyRow = struct(...
        'testCase', "Empty", 'passed', 0, ...
        'model_subspace_id', 0, 'model_id', 0, ...
        'criterionValue', 0, 'pLabels', 0, 'pValues', 0, ...
        'error', "Test not executed, yet.");

for idx = 1:length(names)
    row = emptyRow;
    row.testCase = string(names{idx});
    if idx == 1
        results = row;
    else
        results(idx) = row;
    end
end
results = struct2table(results);

end


% -------------------------------------------------------------------------
function results = evaluateTestResult(results, idx)
% load and compare the selected and expected model .yaml files

% read expectation and results
problem = ReadYaml('petab_select_problem.yaml');
criterion = problem.criterion;
expected = ReadYaml('expected.yaml');
selected = ReadYaml(fullfile('PEtabSelect', 'Results', 'selected_model.yaml'));

% define comparison function
tolerance = 1e-3;
areClose = @(a, b) abs(a - b) < tolerance;

% check IDs
results{idx, "model_subspace_id"} = strcmp(selected.model_subspace_id, expected.model_subspace_id);
results{idx, "model_id"} = strcmp(selected.model_id, expected.model_id);

% check criterion value
results{idx, "criterionValue"} = areClose(selected.criteria.(criterion), expected.criteria.(criterion));

% check the estimated parameters
if isempty(expected.estimated_parameters) == isempty(selected.estimated_parameters)
    if isempty(selected.estimated_parameters)
        results{idx, "pLabels"} = 1;
        results{idx, "pValues"} = 1;
    else
        % sort the fields of the estimated parameters (for easier comparison)
        expected.estimated_parameters = orderfields(expected.estimated_parameters);
        selected.estimated_parameters = orderfields(selected.estimated_parameters);

        % check if the labels of the estimated parameters are the same
        expectedLabels = fieldnames(expected.estimated_parameters);
        selectedLabels = fieldnames(selected.estimated_parameters);
        results{idx, "pLabels"} = isequal(expectedLabels, selectedLabels);

        % check if the values of the estimated parameters are the same
        if results{idx, "pLabels"}
            for i=1:length(expectedLabels)
                results{idx, "pValues"} = areClose( ...
                    selected.estimated_parameters.(selectedLabels{i}), ...
                    expected.estimated_parameters.(expectedLabels{i}));
                if ~results{idx, "pValues"}
                    break
                end
            end
        end
    end
end

% clear error message
results{idx, "error"} = "";
% check if all comparisons are satisfied
results{idx, "passed"} = all(results{1, 3:7}, 'all');
end
