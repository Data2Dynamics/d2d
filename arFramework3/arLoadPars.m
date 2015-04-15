% load model parameters and parameter setting from .mat
% and reconcile with current parameters
%
% arLoadPars(filename, fixAssigned, only_values)
%
% filename:     source file name
%               or number according to the filelist
% fixAssigned:  fix the assigned parameters
% only_values:  only load parameter values, not bound and status
% pfad          path to the results folder, default: './Results'
% 
% 
%  Examples for loading several parameter sets:
% ars = arLoadPars({2,5});
% [ars,ps] = arLoadPars({'Result1','ResultNr2'});
% 
%  Examples for loading all parameter sets:
% ars = arLoadPars('all');
% [ars,ps] = arLoadPars('all');
% arFits(ps)
% 
%  Example:
% arLoadPars('20141112T084549_model_fitted',[],[],'../OtherFolder/Results')

function varargout = arLoadPars(filename, fixAssigned, only_values, pfad)
if(~exist('pfad','var') || isempty(pfad))
    pfad = './Results';
else
    if(strcmp(pfad(end),filesep)==1)
        pfad = pfad(1:end-1);
    end
end

global ar

if(~exist('filename','var') || isempty(filename))
    [~, filename] = fileChooser(pfad, 1, true);
elseif(isnumeric(filename)) % filename is the file-number
    [~, ~, file_list] = fileChooser(pfad, 1, -1);    
    filename = file_list{filename};
elseif(strcmp(filename,'end'))
    filelist = fileList(pfad);
    filename = filelist{end};
elseif(strcmp(filename,'all'))
    filename = fileList(pfad);
%     filename = filename(1:4)
end

if(~exist('fixAssigned', 'var') || isempty(fixAssigned))
    fixAssigned = false;
end
if(~exist('only_values', 'var') || isempty(only_values))
    only_values = false;
end


if(~iscell(filename))    
    ar = arLoadParsCore(ar, filename, fixAssigned, only_values, pfad);
    varargout = cell(0);
else
    ars = cell(size(filename));
    for i=1:length(filename)
        if(isnumeric(filename{i}))
            [~, ~, file_list] = fileChooser(pfad, 1, -1);
            file = file_list{filename{i}};
        else
            file = filename{i};
        end

        ars{i} = arLoadParsCore(ar, file, fixAssigned, only_values, pfad);
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




function ar = arLoadParsCore(ar, filename, fixAssigned, only_values, pfad)
N = 1000;

filename_tmp = filename;

filename = [pfad,'/' filename '/workspace.mat'];
filename_pars = [pfad,'/' filename_tmp '/workspace_pars_only.mat'];
if(exist(filename_pars,'file'))
    S = load(filename_pars);
else
    S = load(filename);
end

fprintf('parameters loaded from file %s:\n', filename);

ass = NaN(size(ar.p));
for j=1:length(ar.p)
    qi = ismember(S.ar.pLabel, ar.pLabel{j});
    
    if(isempty(qi) || sum(qi) == 0)
        ass(j) = 0;
        if(length(ar.p)<=N)
            fprintf('                      %s\n', ar.pLabel{j});
        end
    else
        ass(j) = 1;
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
            if(length(ar.p)<=N)
                fprintf('fixed and assigned -> %s\n', ar.pLabel{j});
            end
        else
            if(length(ar.p)<=N)
                fprintf('          assigned -> %s\n', ar.pLabel{j});
            end
        end
    end
end

if(length(ar.p)>N)
    nnot = length(ass)-sum(ass);
    fprintf('%i assigned, %i not assigned.\n',sum(ass),nnot);
    if(nnot<=30 && nnot>0)
        fprintf('Not assigned are: %s \n',sprintf('%s, ',ar.pLabel{find(ass==0)}));
    end
end

