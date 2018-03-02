% arImportPars(pStruct)
% arImportPars(pStruct, pars_only, pattern, fixAssigned)
% arImportPars(pStruct, pars_only, pattern, fixAssigned, ar)
% 
% 
%   Used by arLoadPars
% 
%   ar      ar-Struct if the parameters should not be imported into the
%           global ar
%           If empty or not provided, then the global ar is used.

function varargout = arImportPars(pStruct, pars_only, pattern, fixAssigned, ar, antipattern)
if(~exist('fixAssigned', 'var') || isempty(fixAssigned))
    fixAssigned = false;
end
if(~exist('pars_only', 'var') || isempty(pars_only))
    pars_only = false;
end
if(~exist('pars_only', 'var') || isempty(pars_only))
    pars_only = false;
end
if(~exist('antipattern', 'var') || isempty(antipattern))
    antipattern = [];
end
if(~exist('pattern', 'var') || isempty(pattern))
    pattern = [];
end
if(~exist('ar', 'var') || isempty(ar))
    global ar
end

N = 1000;  % Number of output message lines.

% remember old values
arIn.p = ar.p;
arIn.lb = ar.lb;
arIn.ub = ar.ub;
arIn.qFit = ar.qFit;
arIn.qLog10 = ar.qLog10;
arIn.type = ar.type;
arIn.mean = ar.mean;
arIn.std = ar.std;

fn = fieldnames(arIn);  % these fields are later checked for changes


if(isempty(pattern))
    js = 1:length(ar.p);
else
    js = find(~cellfun(@isempty,regexp(ar.pLabel, pattern)));
end
if(~isempty(antipattern))
    js(find(~cellfun(@isempty,regexp(ar.pLabel(js), antipattern)))) = [];
end

ass = zeros(size(ar.p));
different = zeros(size(ar.p));
for j=js
    qi = ismember(pStruct.pLabel, ar.pLabel{j});
    
    if(isempty(qi) || sum(qi) == 0)
        ass(j) = 0;
        if(length(ar.p)<=N)
            arFprintf(1, '                      %s\n', ar.pLabel{j});
        end
    else
        ass(j) = 1;
        if(~pars_only)
            ar.p(j) = pStruct.p(qi);
            ar.qLog10(j) = pStruct.qLog10(qi);
            ar.qFit(j) = pStruct.qFit(qi);
            ar.lb(j) = pStruct.lb(qi);
            ar.ub(j) = pStruct.ub(qi);
            if isfield(pStruct,'type')
                ar.type(j) = pStruct.type(qi);
            end
            if isfield(pStruct,'mean')
                ar.mean(j) = pStruct.mean(qi);
            end
            if isfield(pStruct,'std')
                ar.std(j) = pStruct.std(qi);
            end
        else
            if(ar.qLog10(j) == pStruct.qLog10(qi))
                ar.p(j) = pStruct.p(qi);
            elseif(ar.qLog10(j)==1 && pStruct.qLog10(qi)==0)
                ar.p(j) = log10(pStruct.p(qi));
            elseif(ar.qLog10(j)==0 && pStruct.qLog10(qi)==1)
                ar.p(j) = 10^(pStruct.p(qi));
            end
            
            % check bound
            different = (ar.p<ar.lb) | (ar.p>ar.ub);
            ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);
            ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);
        end
        
        newstr = '';
        somethingNew = 0;
        for f=1:length(fn)
            n_char = length(fn{f})+1;
            if abs(ar.(fn{f})(j) - arIn.(fn{f})(j))>1e-10
                if isempty(newstr)
                    newstr = '   new:';
                end
                newstr = sprintf(['%s%',num2str(n_char),'s'],newstr,fn{f});
                somethingNew = 1;
            else
                newstr = sprintf(['%s%',num2str(n_char),'s'],newstr,'');
            end
        end
        if somethingNew == 0
            newstr = '   nothing changed';
        end
        if ( different(j) )
            newstr = sprintf( '%s (clamped against parameter bound)', newstr );
        end
    
        if(fixAssigned)
            ar.qFit(j) = 0;
            if(length(ar.p)<=N)
                arFprintf(1, 'fixed and assigned -> %s%s\n', ar.pLabel{j},newstr);
            end
        else
            if(length(ar.p)<=N)
                arFprintf(1, '          assigned -> %22s%s\n', ar.pLabel{j},newstr);
            end
        end
    end
end


nnot = length(ass)-sum(ass);
if ( nnot > 0 )
    arFprintf(1, '%i parameters were assigned in the destination model (%i not assigned).\n',sum(ass),nnot);
    if(nnot<=30 && nnot>0)
        arFprintf(1, 'Not assigned are: %s \n',sprintf('%s, ',ar.pLabel{ass==0}));
    end
else
    arFprintf(1, 'All parameters assigned.\n');
end

nnot = length(pStruct.p)-sum(ass);
if ( nnot > 0 )
    arFprintf(1, 'There were %i more parameters in the loaded struct than in the target model.\n',nnot);
    tmp=ismember(pStruct.pLabel, ar.pLabel);
    arFprintf(1, 'Not assigned from loaded struct are: %s \n',sprintf('%s, ',pStruct.pLabel{~tmp}));
end

if nargout>0
    varargout{1} = ar;
end