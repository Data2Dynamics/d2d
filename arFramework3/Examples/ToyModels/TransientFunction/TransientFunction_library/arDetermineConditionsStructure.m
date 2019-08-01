% This function tries to determine the conditions in the model:
% 
%   isDefault{m}(ip,c)   n_pold x n_condition 
%                   indicates whether the condition c has default value for
%                   pold{ip}
% 
% 
% Example
%   arLoad % load an ODE model
%   [condsNum,vals,conds,isDefault] = arDetermineConditionsStructure


function [condsNum,vals,conds,isDefault] = arDetermineConditionsStructure
global ar
conds = cell(size(ar.model));
%% initialize conds variable
% this variable will have fielnames pold and contain each fp for all
% conditions
for m=1:length(ar.model)
    conds{m} = struct;
    for c=1:length(ar.model(m).condition)
        for ip=1:length(ar.model(m).condition(c).pold)
            conds{m}.(ar.model(m).condition(c).pold{ip}) = cell(0);
        end
    end
end

for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        for ip=1:length(ar.model(m).condition(c).pold)
            fp = ar.model(m).condition(c).fp(ip);
            fp = regexprep(fp,'^\((\d+)\)$','$1'); % replacement of (0) -> (0) usw.
%             conds{m}.(ar.model(m).condition(c).pold{ip}) = union(conds{m}.(ar.model(m).condition(c).pold{ip}),fp);
            conds{m}.(ar.model(m).condition(c).pold{ip})(end+1) = fp;
        end
    end
end

%% Num find numerival values for the fp stored in conds by performing parameter replacements:
condsNum = conds;

warning off
pval = ar.p;
pval(ar.qLog10==1) = 10.^pval(ar.qLog10==1);
pval = arSym(pval);
pLabel = arSym(ar.pLabel);
for m=1:length(condsNum)
    fn = fieldnames(condsNum{m});        
    for f=1:length(fn)
        indFormula = find(~cellfun(@isempty,regexp(condsNum{m}.(fn{f}),'[a-zA-Z][a-zA-Z]'))); % 2 subsequent letters
        fprintf('Substitutions: Condition-variable %s (%i of %i): %i Formulas found.\n',fn{f},f,length(fn),length(indFormula));
        for i=1:length(indFormula)
            fp = char(arSubs(condsNum{m}.(fn{f}){indFormula(i)},pLabel,pval));
%             fp = regexprep(fp,'^\((\d+)\)$','$1');
            condsNum{m}.(fn{f}){indFormula(i)} = eval(fp);
            if sum(isnan(condsNum{m}.(fn{f}){indFormula(i)}))>0
               condsNum{m}.(fn{f}){indFormula(i)} 
            end
        end        
        notFormula = setdiff(1:length(condsNum{m}.(fn{f})),indFormula);
        for i=1:length(notFormula)
            condsNum{m}.(fn{f}){notFormula(i)} = str2num(condsNum{m}.(fn{f}){notFormula(i)});
            if sum(isnan(condsNum{m}.(fn{f}){notFormula(i)}))>0
               condsNum{m}.(fn{f}){notFormula(i)} 
            end
        end
    end
end
warning on

%% Determine unique values, min, max, default
vals = cell(size(condsNum));
isDefault = cell(size(ar.model));
for m=1:length(condsNum)
    fn = fieldnames(condsNum{m});
    isDefault{m} = zeros(length(fn),length(ar.model(m).condition));
    for f=1:length(fn)
        tmp = condsNum{m}.(fn{f});
        [lev,~,ic] = unique([tmp{:}]);
        anzlev = NaN(1,length(lev));
        for i=1:length(lev)
            anzlev(i) = sum(ic==i);
        end
        
        [~,imax] = max(anzlev);
%         if length(lev)>1
            vals{m}.(fn{f}).default = lev(imax(1));
            vals{m}.(fn{f}).min = min(lev);
            vals{m}.(fn{f}).max = max(lev);
            vals{m}.(fn{f}).anz = length(lev);
            vals{m}.(fn{f}).levels = lev;
            vals{m}.(fn{f}).isDefault = zeros(size(condsNum{m}.(fn{f})));
            % which conditions have default value
            tmp = condsNum{m}.(fn{f});
            isDefault{m}(f,[tmp{:}]==vals{m}.(fn{f}).default) = 1;
%         end
    end
end
