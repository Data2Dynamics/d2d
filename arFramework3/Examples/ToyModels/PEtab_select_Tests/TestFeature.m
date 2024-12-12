% script that runs the petab_select tests

fprintf(2, 'TEST FOR PETAB EXTENSION >>PETAB-SELECT<<\n');

nCases = 9;
tests = 1:nCases;
cases = arrayfun(@(x) sprintf('%04d', x), tests, 'UniformOutput', false);

% results data structures
results = table();
isOk = [];
allExpected = {};
allSelected = {};


global pets
global ar

% set up the paths for the test cases
mainDir = pwd();
testsMainDir = fullfile(mainDir, 'test_cases');
executionMainDir = fullfile(mainDir, 'execution');

%% Run Tests (and evaluate individual results)

try
    parpool();
catch
end

for idx = 1:length(tests)
    % define problem
    testCase = cases{idx};
    testDir = fullfile(testsMainDir, testCase);
    testProblem = fullfile(testDir, 'petab_select_problem.yaml');
    
    % go to execution folder
    executionDir = fullfile(executionMainDir, testCase);
    mkdir(executionDir);
    cd(executionDir);
    
    try
        % run PEtabSelect
        if tests(idx) == 9
            arPetsRunSelect(testProblem, @fitFunction4FinalTest);
        else
            arPetsRunSelect(testProblem, @fitFunction4Tests);
        end

        % load expected and obtained results from files
        selProblem = ReadYaml(testProblem);
        expected = ReadYaml(fullfile(testDir, 'expected.yaml'));
        selected = ReadYaml(fullfile(executionDir, 'selected_model.yaml'));
        allExpected{end+1} = expected;
        allSelected{end+1} = selected;
    
        % evaulate success
        [isOk(end+1), row] = compareModelYamls(expected, selected, selProblem.criterion, testCase);
        results = [results; row];
        if isOk(idx) ~= 1
            fprintf('test case %s failed\n', testCase)
        end

    catch ME
        fprintf('error during test of case %s\n', testCase)
        isOk(end+1) = 2;
        row = struct(...
            'testCase', testCase, 'model_subspace_id', 0, 'model_id', 0, ...
            'criterionValue', 0, 'pLabels', 0, 'pValues', 0, ...
            'error', string(getReport(ME, "extended", "hyperlinks", "off")));
        results = [results; struct2table(row)];
    end
    
end

% after running all tests -> return to main directory
cd(mainDir)

% save the results table
writetable(results, fullfile(mainDir, 'testResults.csv'));

%% Evaluate total results 

if all(isOk==1)
    fprintf(2, 'Correct results all test cases: %s\n', strjoin(cases(isOk==1)));
    fprintf(2, 'PASSED\n' );
else    
    fprintf(2, 'Correct results in test case(s) %s\n', strjoin(cases(isOk==1)));
    fprintf(2, 'Results not identical in test case(s) %s\n', strjoin(cases(isOk==0)));
    fprintf(2, 'Errors in test case(s) %s\n', strjoin(cases(isOk==2)));
    error('FAILED');
end


function [isOk, results] = compareModelYamls(expected, selected, criterion, testCase)
% check model_subspace_id, model_id, criterion and values of estimated parameters
% summarize results in a table and a boolean value

% define comparison function
tolerance = 1e-3;
areClose = @(a, b) abs(a - b) < tolerance;

% initialize the results struct
results = struct();
results.testCase = testCase;
results.model_subspace_id = 0;
results.model_id = 0;
results.criterionValue = 0;
results.pLabels = 0;
results.pValues = 0;
results.error = "";

% check model_subspace_id, model_id
results.model_subspace_id = strcmp(selected.model_subspace_id, expected.model_subspace_id);
results.model_id = strcmp(selected.model_id, expected.model_id);

% check criterion
results.criterionValue = areClose(selected.criteria.(criterion), expected.criteria.(criterion));

% check the estimated parameters
if isempty(expected.estimated_parameters) == isempty(selected.estimated_parameters)
    if isempty(selected.estimated_parameters)
        results.pLabels = 1;
        results.pValues = 1;
    else
        % sort the fields of the estimated parameters (for easier comparison)
        expected.estimated_parameters = orderfields(expected.estimated_parameters);
        selected.estimated_parameters = orderfields(selected.estimated_parameters);

        % check if the labels of the estimated parameters are the same
        expectedLabels = fieldnames(expected.estimated_parameters);
        selectedLabels = fieldnames(selected.estimated_parameters);
        results.pLabels = isequal(expectedLabels, selectedLabels);

        % check if the values of the estimated parameters are the same
        if results.pLabels
            for i=1:length(expectedLabels)
                results.pValues = areClose( ...
                    selected.estimated_parameters.(selectedLabels{i}), ...
                    expected.estimated_parameters.(expectedLabels{i}));
                if ~results.pValues
                    break
                end
            end
        end
    end
end

% create a table with the results
results = struct2table(results);
isOk = all(results{1, 2:6}, 'all');

end


function fitFunction4Tests()

global ar

% deactivate Bessel correction
ar.config.useFitErrorCorrection = 0;

arFit();
arPetsFitSmartly();
arFitLHS(10);

end


function fitFunction4FinalTest()

global ar

% deactivate Bessel correction
ar.config.useFitErrorCorrection = 0;

% use settings similar to pyPESTO
ar.config.rtol = 1e-16;
ar.config.atol = 1e-12;
ar.config.maxsteps = 1e6;

% stricter equilibration settings
ar.config.eq_rtol = 1e-12;
ar.config.eq_tol = 1e-12;
ar.config.init_eq_step = 1e3;

% stricter optimization settings
ar.config.optim.MaxIter = 1e6;
ar.config.optim.TolX = 1e-16;

% change steady state calculation start time (0 instead of -1e7)
arClearEvents();
arSteadyState(1, 1, 1, {}, 0);

% perform fits
arFit();
arPetsFitSmartly();
arFitLHS(40);

% store waterfall plots
close all
arPlotFits();
f = [figure(1), figure(2)];
figurefolder = 'CalibrationPlots';
[~, ~, ~] = mkdir(figurefolder);
savefig(f, fullfile(figurefolder, ar.info.petsModelHash), "compact");
close all

end