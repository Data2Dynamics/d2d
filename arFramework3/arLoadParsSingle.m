% load specific model parameters from .mat
%
% values = arLoadPars(filename, parameternames)
%
% filename:     source file name

function values = arLoadParsSingle(filename, parameternames)

filename_tmp = filename;
filename = ['./Results/' filename '/workspace.mat'];
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
