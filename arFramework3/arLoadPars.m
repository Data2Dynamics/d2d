% load model parameters and parameter setting from .mat
% and reconcile with current parameters
%
% arLoadPars(filename, fixAssigned, only_values)
%
% filename:     source file name
% fixAssigned:  fix the assigned parameters
% only_values:  only load parameter values, not bound and status

function arLoadPars(filename, fixAssigned, only_values)

global ar

if(~exist('filename','var') || isempty(filename))
    [~, filename] = fileChooser('./Results', 1, true);
else
    if(strcmp(filename,'end'))
        filelist = fileList('./Results');
        filename = filelist{end};
    end
end
filename_tmp = filename;
filename = ['./Results/' filename '/workspace.mat'];
filename_pars = ['./Results/' filename_tmp '/workspace_pars_only.mat'];

if(~exist('fixAssigned', 'var'))
    fixAssigned = false;  
end
if(~exist('only_values', 'var'))
    only_values = false;  
end

if(exist(filename_pars,'file'))
    S = load(filename_pars);
    fprintf('parameters loaded from file %s:\n', filename_pars);
else
    S = load(filename);
    fprintf('parameters loaded from file %s:\n', filename);
end

for j=1:length(ar.p)
    qi = ismember(S.ar.pLabel, ar.pLabel{j});
    
    if(isempty(qi) || sum(qi) == 0)
        fprintf('                      %s\n', ar.pLabel{j});
    else
        if(~only_values)
            ar.p(j) = S.ar.p(qi);
            ar.qLog10(j) = S.ar.qLog10(qi);
            ar.qFit(j) = S.ar.qFit(qi);
            ar.lb(j) = S.ar.lb(qi);
            ar.ub(j) = S.ar.ub(qi);
            if isfield(S.ar,'type')
                ar.type(j) = S.ar.type(qi);
            end
            if isfield(S.ar,'mean')
                ar.mean(j) = S.ar.mean(qi);
            end
            if isfield(S.ar,'std')
                ar.std(j) = S.ar.std(qi);
            end
        else
            if(ar.qLog10(j) == S.ar.qLog10(qi))
                ar.p(j) = S.ar.p(qi);
            elseif(ar.qLog10(j)==1 && S.ar.qLog10(qi)==0)
                ar.p(j) = log10(S.ar.p(qi));
            elseif(ar.qLog10(j)==0 && S.ar.qLog10(qi)==1)
                ar.p(j) = 10^(S.ar.p(qi));
            end
            
            % check bound
            ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);
            ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
        end
        
        if(fixAssigned)
            ar.qFit(j) = 0;
            fprintf('fixed and assigned -> %s\n', ar.pLabel{j});
        else
            fprintf('          assigned -> %s\n', ar.pLabel{j});
        end
    end
end
