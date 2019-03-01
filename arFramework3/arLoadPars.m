% [ars,ps] = arLoadPars([resultFolder], [fixAssigned], [pars_only], [result_path], [pattern], [antipattern])
%
% load model parameters and parameter setting from
% result_path/resultFolder/workspace_pars_only.mat and reconcile with
% current parameters 
%
%   resultFolder  result folder name(s)                       []
%                 or number according to the filelist
%                 or a previously created global variable ar
%                 or uses fileChooser to allow user selection 
%   fixAssigned   fix the assigned parameters                 [false]
%   pars_only     only load parameter values, not bound and status [false]
%   result_path   path to the results folder                  ['./Results']
%   pattern       search pattern within parameter names       []
%   antipattern   pattern for parameters to exclude from load []
% 
%   ars           cell array of workspaces with loaded parameters (optional output)
%   ps            array of parameter sets (optional output)
%  
% Examples 
%   loading several parameter sets:
%   >> ars = arLoadPars({2,5});
%   >> [ars,ps] = arLoadPars({'Result1','ResultNr2'});
% 
%   loading all parameter sets:
%   >> ars = arLoadPars('all');
%   >> [ars,ps] = arLoadPars('all');
%   >> arFits(ps)
% 
%   loading from different source folder:
%   >> arLoadPars('20141112T084549_model_fitted',[],[],'../OtherFolder/Results')
% 
%   loading from another ar struct:
%   >> tmp = load(['testfile.mat']);  % testfile.mat has to contain an ar struct
%   >> arLoadPars(tmp.ar)

% 
% See also arLoadParsSingle, arImportPars

function varargout = arLoadPars(resultFolder, fixAssigned, pars_only, result_path, pattern, antipattern)
if(~exist('result_path','var') || isempty(result_path))
    result_path = './Results';
else
    if(strcmp(result_path(end),filesep)==1)
        result_path = result_path(1:end-1);
    end
end

global ar

if(~exist('fixAssigned', 'var') || isempty(fixAssigned))
    fixAssigned = false;
end
if(~exist('pars_only', 'var') || isempty(pars_only))
    pars_only = false;
end
if(~exist('pattern', 'var') || isempty(pattern))
    pattern = [];
end
if(~exist('antipattern', 'var' ) || isempty(antipattern))
    antipattern = [];
end

if(~exist('resultFolder','var') || isempty(resultFolder))
    [~, resultFolder] = fileChooser(result_path, 1, true);
elseif(isnumeric(resultFolder)) % resultFolder is the file-number
    [~, ~, file_list] = fileChooser(result_path, 1, -1);    
    resultFolder = file_list{resultFolder};
elseif(strcmp(resultFolder,'end'))
    filelist = fileList(result_path);
    resultFolder = filelist{end};
elseif(strcmp(resultFolder,'all'))
    resultFolder = fileList(result_path);
elseif ischar(resultFolder)
    if ~exist(result_path,'dir')
        error('Folder %s does not exist.',result_path);
    end
    [~,resultFolder]=fileparts(resultFolder);    % remove path
end


if(~iscell(resultFolder))    
    ar = arLoadParsCore(ar, resultFolder, fixAssigned, pars_only, result_path, pattern, antipattern);
    varargout = cell(0);
else
    ars = cell(size(resultFolder));
    for i=1:length(resultFolder)
        if(isnumeric(resultFolder{i}))
            [~, ~, file_list] = fileChooser(result_path, 1, -1);
            file = file_list{resultFolder{i}};
        else
            file = resultFolder{i};
        end

        ars{i} = arLoadParsCore(ar, file, fixAssigned, pars_only, result_path, pattern, antipattern);
    end
    varargout{1} = ars;
    if nargout>1
        ps = NaN(length(ars),length(ar.p));
        for i=1:length(ars)
            ps(i,:) = ars{i}.p;
        end
        varargout{2} = ps;
    end
end

% Invalidate cache so simulations do not get skipped
arCheckCache(1);

function ar = arLoadParsCore(ar, resultFolder, fixAssigned, pars_only, pfad, pattern, antipattern)

if(ischar(resultFolder))    
    filename_pars = [pfad,'/' resultFolder '/workspace.mat'];
    filename_pars_only = [pfad,'/' resultFolder '/workspace_pars_only.mat'];
    
    if(exist(filename_pars_only,'file'))
        S = load(filename_pars_only);
    else
        S = load(filename_pars);
    end
    arFprintf(1, 'parameters loaded from file %s:\n', resultFolder);
elseif(isstruct(resultFolder))
    S = struct([]);
    S(1).ar = resultFolder;
else
    error('not supported variable type for resultFolder');
end

ar = arImportPars(S.ar, pars_only, pattern, fixAssigned, ar, antipattern);
