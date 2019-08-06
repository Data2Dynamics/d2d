function y = arExportPEtab(dummy)
% function y = arExportPEtab
%
% Exports data, condition and measurement details to PEtab format
% See also https://github.com/ICB-DCM/PEtab

% Not implemented yet:
% - Measurement table: preEquilibrationConditionId,
% - Visualization table
% - Parameter table: priorType, priorParameters: not yet specified!

global ar

%% Condition and Measurment Table
% Collect conditions
condPar = {}; condVal = {}; conditionID = {};
for imodel = 1:length(ar.model)
    for idata = 1:length(ar.model(imodel).data)
        for icond = 1:length(ar.model(imodel).data(idata).condition)
            condPar{end+1} = ar.model(imodel).data(idata).condition(icond).parameter;
        end
    end
end
peConds = unique(condPar);

% Write condition and measurement tables
measT = table;
peCondValues = [];
numOfConds = length(peConds);
conditionID = {};

for imodel = 1:length(ar.model)
    for idata = 1:length(ar.model(imodel).data)
        
        % conditions
        condPar = {}; condVal = []; condPos = []; conditionID_tmp = '';
        for icond = 1:length(ar.model(imodel).data(idata).condition)
            
            condPar{end+1} = ar.model(imodel).data(idata).condition(icond).parameter;
            condVal = [condVal str2num(ar.model(imodel).data(idata).condition(icond).value)];
            condPos = [condPos find(contains(peConds, condPar{icond}))];
            
            %conditionID_tmp = [conditionID_tmp '_' ar.model(imodel).data(idata).condition(icond).parameter ...
            %'_' ar.model(imodel).data(idata).condition(icond).value];
        end
        
        rowToAdd = zeros(1, numOfConds);
        rowToAdd(condPos) = condVal;
        peCondValues(end+1,:) = rowToAdd;
        %conditionID{end+1} = conditionID_tmp;
        simuConditionID = ['model' num2str(imodel) '_data' num2str(idata)];
        conditionID{end+1} = simuConditionID;
        
        % write data
%         observableIDs = repmat(ar.model(imodel).data(idata).y, ...
%             [length(ar.model(imodel).data(idata).yExp) 1]);
%         observableIDs = observableIDs(:);
%                 
%         measurements = ar.model(imodel).data(idata).yExp(:);
%           
%         timepoints = repmat(ar.model(imodel).data(idata).tExp, [size(ar.model(imodel).data(idata).yExp, 2) 1]);
% 
%         simulationConditionIDs = repmat(simuConditionID, [length(timepoints) 1]);
% 
%         obsPars_tmp = cellfun(@strsplit, ar.model(imodel).data(idata).fy, repmat({{'+','-','*','/','(',')','^',' '}}, [length(ar.model(imodel).data(idata).fy) 1]), 'UniformOutput', false);
%         obsPars_tmp = cellfun(@intersect, obsPars_tmp, repmat({ar.pLabel},[length(ar.model(imodel).data(idata).fy) 1]), repmat({'stable'},[length(ar.model(imodel).data(idata).fy) 1]), 'UniformOutput', false);
%         observableParameters = obsPars_tmp;
%         
%         noisePars_tmp = cellfun(@strsplit, ar.model(imodel).data(idata).fystd', repmat({{'+','-','*','/','(',')','^',' '}}, [length(ar.model(imodel).data(idata).fystd) 1]), 'UniformOutput', false);
%         noisePars_tmp = cellfun(@intersect, noisePars_tmp, repmat({ar.pLabel},[length(ar.model(imodel).data(idata).fystd) 1]), repmat({'stable'},[length(ar.model(imodel).data(idata).fystd) 1]), 'UniformOutput', false);
%         noiseParameters = noisePars_tmp;
        
        % measurements
        for iy = 1:length(ar.model(imodel).data(idata).y)
            observableIDs = repmat(ar.model(imodel).data(idata).y(iy), ...
                [length(ar.model(imodel).data(idata).yExp(:,iy)) 1]);
            
            measurements = ar.model(imodel).data(idata).yExp(:, iy);
            timepoints = ar.model(imodel).data(idata).tExp;
            simulationConditionIDs = repmat({simuConditionID}, [length(timepoints) 1]);
            
            obsPars_tmp = strsplit(ar.model(imodel).data(idata).fy{iy}, {'+','-','*','/','(',')','^',' '});
            obsPars_tmp = intersect(obsPars_tmp, ar.pLabel);
            obsPars_tmp = strcat(obsPars_tmp, ';');
            observableParameters = repmat({[obsPars_tmp{:}]}, [length(timepoints) 1]);
            
            noisePars_tmp = strsplit(ar.model(imodel).data(idata).fystd{iy}, {'+','-','*','/','(',')','^',' '});
            noisePars_tmp = intersect(noisePars_tmp, ar.pLabel);
            noisePars_tmp = strcat(noisePars_tmp, repmat(';', [length(noisePars_tmp)-1 1]));
            noiseParameters = repmat({[noisePars_tmp{:}]}, [length(timepoints) 1]);
            
            obsTrafos = cell(length(timepoints),1);
            obsTrafos(:) = {'lin'}; obsTrafos(logical(ar.model(imodel).data(idata).logfitting(iy))) = {'log'};
            
            noiseDist = cell(length(timepoints),1);
            noiseDist(:) = {'normal'};
            
            measT = [measT; table(observableIDs, measurements, timepoints, ...
                simulationConditionIDs, observableParameters, noiseParameters, obsTrafos, noiseDist)];
        end
    end
end

condT = array2table(peCondValues);
condT.Properties.VariableNames = peConds;

condT = [table(conditionID') condT];
condT.Properties.VariableNames{1} = 'conditionId';

measT.Properties.VariableNames = {'observableId', 'simulationConditionId', ...
    'measurement', 'time', 'observableParameters', 'noiseParameters',...
    'observableTransformation', 'noiseDistribution'};

writetable(condT, 'peTAB_conditions.tsv', 'Delimiter', '\t', 'FileType', 'text')
writetable(measT, 'peTAB_measurments.tsv', 'Delimiter', '\t', 'FileType', 'text')

%% Parameter Table
parameterScale_tmp = cell(1, length(ar.qLog10));
parameterScale_tmp(:) = {'lin'};
parameterScale_tmp(ar.qLog10 == 1) = {'log10'};

nominalValue_tmp = ar.p;
%nominalValue_tmp = ar.qLog10.*nominalValue_tmp + (1-ar.qLog10).*10.^nominalValue_tmp;

lowerBound_tmp = ar.lb;
%lowerBound_tmp = ar.qLog10.*lowerBound_tmp + (1-ar.qLog10).*10.^lowerBound_tmp;

upperBound_tmp = ar.ub;
%upperBound_tmp = ar.qLog10.*upperBound_tmp + (1-ar.qLog10).*10.^upperBound_tmp;

estimate_tmp = ar.qFit;
estimate_tmp(ar.qFit == 2) = 0;

parameterID = ar.pLabel;
parameterName = ar.pLabel;
parameterScale = parameterScale_tmp;
lowerBound = lowerBound_tmp;
upperBound = upperBound_tmp;
nominalValue = nominalValue_tmp;
estimate = estimate_tmp;

parT = table(parameterID(:), parameterName(:), parameterScale(:), ...
    lowerBound(:), upperBound(:), nominalValue(:), estimate(:));

parT.Properties.VariableNames = {'parameterID', 'parameterName', ...
    'parameterScale', 'lowerBound', 'upperBound', 'nominalValue', 'estimate',};

writetable(parT, 'peTAB_parameters.tsv', 'Delimiter', '\t', 'FileType', 'text')

end