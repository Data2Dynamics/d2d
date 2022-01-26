fprintf( 2, 'TEST FOR PETAB EXTENSION >>PETAB-SELECT<<\n' );
clear isOk expected actual selProblem Ncases

% do all
cases = {'0001', '0002', '0003', '0004', '0005', '0006', '0007', '0008'};
Ncases = numel(cases);

isOk = NaN(Ncases,1);

%% Check the python environment of system command 
venvActPath = [];
syscom = [venvActPath, 'petab_select --help'];
[status1,~] = system(syscom);
if status1 ~= 0
    venvActPath = '~/d2d_python_venv/bin/activate';
    initstr = sprintf('source %s; ', venvActPath);

    syscom = [initstr, ' petab_select --help'];
    [status2,~] = system(syscom);
    if status2 ~= 0
        error(sprintf('Calling petab_select from the command line failed. ...\nPlease check your Python environment and the PEtab-select installation.'))
    end
end

%parpool(4)

for i = 1:Ncases
    cd(cases{i})
    
    [status,msg,messageid] = rmdir('output','s');
    if status==0
        if  ~strcmp('MATLAB:RMDIR:NotADirectory',messageid)
            fprintf( 2, 'Error while removing folder output.\n Error message: %s\n', msg);
        end
    end
   
    try
        arPEtabSelect(venvActPath)
        
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
    fprintf( 2, 'Correct results all test cases: %s\n', strjoin(cases(logical(isOk(isOk==1)))));
    fprintf( 2, 'PASSED\n' );
else
    isOk(isnan(isOk)) = 2;
    
    fprintf( 2, 'Correct results in test case(s) %s\n', strjoin(cases(logical(isOk==1))));
    fprintf( 2, 'Results not identical in test case(s) %s\n', strjoin(cases(logical(isOk==0))));
    fprintf( 2, 'Errors in test case(s) %s\n', strjoin(cases(logical(isOk==2))));
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

if ~isempty(expected.estimated_parameters) 
    if ~isempty(actual.estimated_parameters)
        names = fieldnames(expected.estimated_parameters);
        for ipar = 1:length(expected.estimated_parameters)
            if round(actual.estimated_parameters.(names{ipar}), 3) ~= round(expected.estimated_parameters.(names{ipar}), 3)
                isOk = 0;
            end
        end
    end
end
end