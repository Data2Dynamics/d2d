fprintf( 2, 'TEST FOR PETAB EXTENSION >>PETAB-SELECT<< M\n' );
clear isOk expected actual selProblem Ncases

% do all
cases = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008'};
Ncases = numel(cases);

isOk = NaN(Ncases,1);
%parpool(4)
for i = 1:Ncases
    cd(cases{i})
    system('rm -r output')
    try
    arPEtabSelect('~/_d2d_python_venv/bin/activate')    
    
    expected = ReadYaml('expected.yaml');
    actual = ReadYaml('petab-select/selected_model.yaml');
    selProblem = ReadYaml('petab_select_problem.yaml');
    isOk(i) = comparePeTabYamlStruct(expected, actual, selProblem.criterion);
    if isOk(i) ~= 1
        fprintf('test case %s failed\n', cases{i})
    end
    catch
        fprintf('error during test of case %s\n', cases{i})
    end
    cd ..
end

if sum(isOk) == Ncases
    fprintf( 2, 'PASSED\n' );
else
    fprintf( 2, 'Errors in test case(s) %s\n', strjoin(cases(logical(~Working)),', '));
    error( 'FAILED');
end

function isOk = comparePeTabYamlStruct(expected, actual, criterion)
% check criteria, model_subspace_id and values of estimated parameters

isOk = 1;
if round(actual.criteria.(criterion),3) ~= round(expected.criteria.(criterion), 3)
    isOk = 0;
end
if actual.model_subspace_id ~= expected.model_subspace_id
    isOk = 0;
end

names = fieldnames(expected.estimated_parameters);
for ipar = 1:length(expected.estimated_parameters)
    if round(actual.estimated_parameters.(names{ipar}), 3) ~= round(expected.estimated_parameters.(names{ipar}), 3)
        isOk = 0;
    end
end
end