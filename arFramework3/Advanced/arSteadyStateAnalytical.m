function varargout = arSteadyStateAnalytical(replacements, ODESSargs, pythonsymlink)

% [success, steadystate, fullout] = arSteadyStateAnalytical(ptyhonsymlink, ODESSargs, replacements)
%
% Calculates analytical steady state from compiled models using the python
% tool ODESS [1,2]. Results are printed to the console and must be manually
% copied into the CONDITIONS section of the model definition file. 
%
%   replacements    Cell array of string of parameters that will be 
%                   replaced with zeros before calculation of the steady
%                   state. May be used to switch off certain inputs/stimuli 
%   ODESSargs       Cell array of arguments passed to ODESS. Prodivide the 
%                   first 3 arguments as strings (e.g. givenCQs = 
%                   "'AKT+pATK=totalAKT', 'MEK+pMEK=totalMEK'"), the fourth as numeric. 
%                   {injections,givenCQs,neglect,sparsifyLevel}
%                   [{'','','',2}]
%   pythonsymlink   Symlink for calling Python 3 from the terminal ['python3']
% 
%
%   success         Boolean variable indicating whether steady state was
%                   found
%   steadystate     Cell array of steady state expressions in the D2D
%                   syntax
%   fullout         Full output of ODESS
%
% Requirements:
%   - Python 3 installation
%   - Python modules sympy, numpy, csv, random 
%
% Steady states can alternatively be obtained numerically by
% pre-equilibrating the model using arSteadyState
%
% References:
% [1] M. Rosenblatt, J. Timmer, D. Kaschek. 
% Customized steady-state constraints for parameter estimation in non-
% linear ordinary differential equation models. 
% Frontiers in Cell and Developmental Biology 4, 2016, 41
% [2] http://sysbio.uni-freiburg.de/mrosen/
% 
% See also: arSteadyState

global ar

if(~exist('replacements','var') || isempty(replacements))
    replacements = '';
end
if(~exist('pythonsymlink','var') || isempty(pythonsymlink))
    pythonsymlink = '';
end
if(~exist('ODESSargs','var') || isempty(ODESSargs))
    ODESSargs = {'','','',2};
end


% check python version
pythonsymlink = 'python3';
[status, cmdout] = system([pythonsymlink ' --version']);
pythonvers = regexp(cmdout, '\d*', 'Match');
if str2num(pythonvers{3}) < 3
    error('Python 3 required for ODESS (analytical calculation of steady states)')
end

for imodel = 1:length(ar.model)
    % write stochiometric matrix
    fprintf('Writing stochiometric matrix for model %i to %s...\n', imodel, [ar.model(imodel).name, '__model.csv'])
    arExportModelToDmod('model',imodel,[],[],replacements);

    % copy ODESS
    copyfile([ar.info.ar_path, filesep, 'ThirdParty', filesep, 'ODESS.py'], pwd)
    
    % write python script
    workFile = fopen('doWork_ODESS.py', 'w');
    fprintf(workFile, 'from ODESS import *\n');
    fprintf(workFile, 'ODESS("%s",[%s],[%s],[%s],%i,"M")', [ar.model(imodel).name, '__model.csv'],...
        ODESSargs{1}, ODESSargs{2}, ODESSargs{3}, ODESSargs{4});
    fclose(workFile);
    
    % do python work
    fprintf('Calculating steady state for model %i...\n', imodel)
    [status, fullout] = system(sprintf('%s %s', pythonsymlink, 'doWork_ODESS.py'));
    try
        steadystate = regexp(fullout, 'Testing Steady State...', 'split');
        steadystate = regexp(steadystate{2}, 'I obtained the following equations:', 'split');
        
        if ~(strtrim(steadystate{1}) == 'Solution is correct!')
            warning('ODESS failed for model %i.\n\n', imodel)
            fprintf('ODESS output:\n%s\n', fullout)
            success = 0;
        else
            fprintf('ODESS successfull for model %i.\n', imodel)
            success = 1;
        end
        
        steadystate = regexp(steadystate{2}, 'Number of Species', 'split');
        steadystate = strtrim(steadystate{1});
        steadystate = regexprep(steadystate,'\n\n','\n');
        steadystate = regexprep(steadystate,'\t','');
        
        fprintf('\n')
        fprintf('ANALYTICAL STEADY STATE:\n')
        fprintf('(Copy into CONDITIONS section of model-def file)\n')
        fprintf('#####################\n')
        fprintf('%s\n', steadystate)
        fprintf('#####################\n')

        
    catch
        warning('ODESS failed for model %i.\n\n', imodel)
        fprintf('ODESS output:\n%s\n', fullout)
        fprintf(2,'Could not calculate analytical steady state for model %i.\n\n', imodel)
        success = 0;
    end
    varargout{1} = success;
    varargout{2} = steadystate;
    varargout{3} = fullout; 
end
end

