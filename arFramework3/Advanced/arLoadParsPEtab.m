%% arLoadParsPEtab(filename)
%
% Loads parameters of parameters_*.tsv into ar
%
% arLoadParsPEtab('parameters_Boehm_JProteomeRes2014.tsv')

function arLoadParsPEtab(filename)

global ar

T = tdfread(filename);

%% Initialize all, in case parameters are not set in parameters.tsv
ar.p = ar.pExtern;
ar.pLabel = ar.pExternLabels;
ar.lb = ones(size(ar.p))*(-5);
ar.ub = ones(size(ar.p))*3;
ar.qLog10 = zeros(size(ar.p));
ar.qFit = zeros(size(ar.p));
ar.qDynamic = ones(size(ar.p));
ar.qError = zeros(size(ar.p));
ar.qInitial = zeros(size(ar.p));
ar.type = zeros(size(ar.p));

%% Get indices of existing parameters to replace
ParsAll = ar.pExternLabels;
ParsRead = cellstr(T.parameterId);
if any(contains(ParsRead,'sd'))
    ParsRead = strrep(ParsRead,'sd','noiseParameter1');  %% OBservable names still false
end
idx=[];
for i=1:size(ParsRead,1)
    if any(contains(ParsAll,ParsRead(i,:)))
        [~,idx(end+1)] = find(ismember(ParsAll,ParsRead(i,:)));
    end
end

%% Set Parameters to info in parameter*.tsv
ar.p(idx) = T.nominalValue;
ar.lb(idx) = T.lowerBound;
ar.ub(idx) = T.upperBound;
ar.qFit(idx) = T.estimate;

for i=1:size(T.parameterScale,1)
    % Label
    ar.pLabel(idx(i)) = {strtrim(T.parameterId(i,:))};
    % qLog10
    if contains(T.parameterScale(i,:),'lin')
        ar.qLog10(i) = 0;
    elseif strcmp(T.parameterScale(i,:),'log10')
        ar.qLog10(i) = 1;
    else
        warning([T.parameterScale(i,:) ' is not yet implemented as qLog. ar.qLog10 is set to 0.'])
        ar.qLog10(i) = 0;
    end
    % Prior
    if ischar(T.priorType(i))
        if isnumeric(T.priorParameters)
            arSetPrior(i,T.priorType(i),T.priorParameters(i,1),T.priorParameters(i,2))
        else
            arSetPrior(i,T.priorType(i))
        end
    end
end

