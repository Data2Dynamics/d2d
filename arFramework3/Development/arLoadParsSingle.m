% values = arLoadParsSingle(resultFolder, parameternames)
%
% get specific model parameters from pwd/Results/resultFolder/workspace_pars_only.mat
% of (if not existing) from pwd/Results/resultFolder/workspace.mat
%
% resultFolder      source file name. This is a folder in the Results
%                   folder
% parameternames    cell array of parameternames that are searched
% 
% values            array of parameter values corresponding.
% 
% Note: The ar struct is not affected by this function call.
%
% see also arLoadPars

function values = arLoadParsSingle(resultFolder, parameternames)
if ~exist('Results','dir')
    error('No results folder exist. arLoadParsSingle can only be executed in a D2D working directory.')
end

filename_tmp = resultFolder;
filename = ['./Results/' resultFolder '/workspace.mat'];
filename_pars = ['./Results/' filename_tmp '/workspace_pars_only.mat'];

if(exist(filename_pars,'file'))
    S = load(filename_pars);
    fprintf('parameters loaded from file %s:\n', filename_pars);
else
    S = load(filename);
    fprintf('parameters loaded from file %s:\n', filename);
end

values = nan(size(parameternames));
for j=1:length(parameternames)
    qi = ismember(S.ar.pLabel, parameternames{j});
    
    if(~isempty(qi) && sum(qi) == 1)

            values(j) = S.ar.p(qi);
            
% TODO implement other information such as:
%             ar.qLog10(j) = S.ar.qLog10(qi);
%             ar.qFit(j) = S.ar.qFit(qi);
%             ar.lb(j) = S.ar.lb(qi);
%             ar.ub(j) = S.ar.ub(qi);
%             if isfield(S.ar,'type')
%                 ar.type(j) = S.ar.type(qi);
%             end
%             if isfield(S.ar,'mean')
%                 ar.mean(j) = S.ar.mean(qi);
%             end
%             if isfield(S.ar,'std')
%                 ar.std(j) = S.ar.std(qi);
%             end
    end
end
