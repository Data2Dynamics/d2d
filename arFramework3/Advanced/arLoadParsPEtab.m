% arLoadParsPEtab(filename)
%
%   filename    name of the PEtab parameter table
%
% See also
%   arImportPEtab

function arLoadParsPEtab(filename)
global ar

if ~contains(filename,'.tsv')
    if ~contains(filename,'.')
        filename = [filename '.tsv'];
    else
        error('this file type is not supported!')
    end
end

T = tdfread(filename);
try
    T.parameterId = cellstr(T.parameterId);
    [BothPars, ia, ib] = intersect(ar.pLabel,T.parameterId);
catch
    T.parameterID = cellstr(T.parameterID);
    [BothPars, ia, ib] = intersect(ar.pLabel,T.parameterID);
end

T.qLog10 = nan(size(T.nominalValue));
T.qLog10(contains(string(T.parameterScale),'lin')) = 0;
T.qLog10(contains(string(T.parameterScale),'log10')) = 1;

if any(isnan(T.qLog10))
    warning([T.parameterScale(isnan(T.parameterScale)) ' is not yet implemented as qLog. ar.qLog10 is set to 0.'])
end

% apply log10 trafo if flag is set. (Values given are on lin scale)
for i = 1:length(ib)
    if T.qLog10(ib(i))
        % arSetPars(pLabel, [p], [qFit], [qLog10], [lb], [ub], [type], [meanp], [stdp])
        arSetPars(ar.pLabel(ia(i)), log10(T.nominalValue(ib(i))), T.estimate(ib(i)), T.qLog10(ib(i)),log10(T.lowerBound(ib(i))),log10( T.upperBound(ib(i))))
    else
        arSetPars(ar.pLabel(ia(i)), T.nominalValue(ib(i)), T.estimate(ib(i)), T.qLog10(ib(i)),T.lowerBound(ib(i)), T.upperBound(ib(i)))
    end
end


% Priors
if isfield(T,'objectivePriorType')

    % Remap PEtab prior types to d2d prior types
    T.objectivePriorType = cellstr(T.objectivePriorType);  % Convert char array to cell array of strings
    petab_to_d2d = containers.Map({'parameterScaleUniform', 'normal', 'uniform', 'laplace'}, {0, 1, 2, 3});
    T.objectivePriorType = cellfun(@(x) petab_to_d2d(x), T.objectivePriorType);


    for i=1:length(BothPars)
        if ~isempty(T.objectivePriorType(ib(i)))
            if T.objectivePriorType(ib(i)) ~= 0
                PriorParameters = str2num(T.objectivePriorParameters(ib(i),:));
                arSetPrior(ia(i),T.objectivePriorType(ib(i),:), ...
                    PriorParameters(1), PriorParameters(2))
            else
                arSetPrior(ia(i),T.objectivePriorType(ib(i),:))
            end
        end
    end
end