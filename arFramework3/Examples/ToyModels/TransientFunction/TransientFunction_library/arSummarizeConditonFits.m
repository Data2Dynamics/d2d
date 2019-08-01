% This function tries to determine the conditions in the model:
% - Which parameters are changed over conditions?
% - How many levels are there?

function [conds,out] = arDetermineConditionsStructure(fits)
global ar
%% How many values have pold?
anz = struct;
for i=1:length(fits)
    for ip=1:length(fits{i}.condition.pold)
        anz.(fits{i}.condition.pold{ip}) = cell(0);
    end
end

for i=1:length(fits)
    for ip=1:length(fits{i}.condition.pold)
        anz.(fits{i}.condition.pold{ip}) = union(anz.(fits{i}.condition.pold{ip}),fits{i}.condition.fp{ip});
    end
end

%%
fn = fieldnames(anz);
severals = fn(cellfun(@length,struct2cell(anz))>1);

out = struct;
out.m = NaN(length(fits),1);
out.x = cell(length(fits),1);
out.c = NaN(length(fits),1);
for s=1:length(severals)
    out.(severals{s}) = cell(length(fits),1);
end

for i=1:length(fits)
    out.m(i) = fits{i}.m;
    out.c(i) = fits{i}.c;
    out.x{i,1} = fits{i}.x;
    for ip=1:length(fits{i}.pLabel)
        out.(fits{i}.pLabel{ip})(i,1) = fits{i}.p(ip);
    end
    out.toffset_TF2(i,1) = out.toffset_TF(i)/10*max(fits{1}.tFine);
    [~,ia,ib] = intersect(severals,fits{i}.condition.pold);
    for ii=1:length(ia)
        out.(severals{ia(ii)}){i,1} = fits{i}.condition.fp{ib(ii)};
    end
end

warning off
pval = ar.p;
pval(ar.qLog10==1) = 10.^pval(ar.qLog10==1);
pval = arSym(pval);
pLabel = arSym(ar.pLabel);
for s=1:length(severals)
    indFormula = find(~cellfun(@isempty,regexp(out.(severals{s}),'[a-zA-Z][a-zA-Z]'))); % 2 subsequent letters
    fprintf('Substitutions: Condition-variable %s (%i of %i): %i Formulas found.\n',severals{s},s,length(severals),length(indFormula));
    for i=1:length(indFormula)
        out.(severals{s}){indFormula(i)} = char(arSubs(out.(severals{s}){indFormula(i)},pLabel,pval));
    end
end
warning on
WriteTable('ConditionStructure.xls',out,[],[],',');



%%
fn = severals;
conds = struct;
for f=1:length(fn)
%     [lev,~,ic] = unique(out.(fn{f}));
    tmp = cellfun(@str2num,out.(fn{f}),'UniformOutput',false);
    [lev,~,ic] = unique([tmp{:}]);
    anzlev = NaN(1,length(lev));
    for i=1:length(lev)
        anzlev(i) = sum(ic==i);
    end
    
    [~,imax] = max(anzlev);
    if length(lev)>1
        conds.(fn{f}).default = lev(imax(1));
        conds.(fn{f}).min = min(lev);
        conds.(fn{f}).max = max(lev);
        conds.(fn{f}).anz = length(lev);
        conds.(fn{f}).levels = lev;
    end
end

