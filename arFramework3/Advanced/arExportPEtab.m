function y = arExportPEtab(dummy)
% function y = arExportPEtab
%
% Export data, condition and measurement details to PEtab data format
% See also https://github.com/ICB-DCM/PEtab

% Not implemented yet:
% - Measurement table: preEquilibrationConditionId,
% - Visualization table
% - Parameter table: priorType, priorParameters: not yet specified!
% - Different noise distributions

global ar

%% Write Export Directory
if(~exist('./SBML_Export', 'dir'))
    mkdir('./SBML_Export')
end

for imodel = 1:length(ar.model)
    %% Condition and Measurement Table
    % Collect conditions
    condPar = {}; condVal = {}; conditionID = {};
    for idata = 1:length(ar.model(imodel).data)
        for icond = 1:length(ar.model(imodel).data(idata).condition)
            condPar{end+1} = ar.model(imodel).data(idata).condition(icond).parameter;
        end
        % catch d2d fix for predictors that are not time
        if ~(strcmp(ar.model(imodel).data(idata).t, 'time') || strcmp(ar.model(imodel).data(idata).t, 't'))
            for it = 1:length(ar.model(imodel).data(idata).tExp)
                condPar{end+1} = [ar.model(imodel).data(idata).t...
                    '_' num2str(it)];
            end
        end
    end
    peConds = unique(condPar);
    
    % Write condition and measurement tables
    measT = table;
    peCondValues = [];
    numOfConds = length(peConds);
    conditionID = {};
    
    for idata = 1:length(ar.model(imodel).data)
        % conditions
        condPar = {}; condVal = []; condPos = []; conditionID_tmp = '';
        
%         % catch d2d fix for predictors that are not time
%         if ~(strcmp(ar.model(imodel).data(idata).t, 'time') || strcmp(ar.model(imodel).data(idata).t, 't'))
%             for it = 1:length(ar.model(imodel).data(idata).tExp)
%                 condPar{end+1} = [ar.model(imodel).data(idata).t...
%                     '_' num2str(it)];
%                 condVal = [condVal ar.model(imodel).data(idata).tExp(it)];
%                 condPos = [condPos find(strcmp(peConds, condPar{it}))];
%             end
%         end
        
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
                
        % measurements
        for iy = 1:length(ar.model(imodel).data(idata).y)
            observableIDs = repmat(ar.model(imodel).data(idata).y(iy), ...
                [length(ar.model(imodel).data(idata).yExp(:,iy)) 1]);
            
            measurements = ar.model(imodel).data(idata).yExp(:, iy);
            timepoints = ar.model(imodel).data(idata).tExp;
            
            obsTrafos = cell(length(timepoints),1);
            obsTrafos(:) = {'log10'};
            
            if ~ar.model(imodel).data(idata).logfitting(iy)
                obsTrafos(:) = {'lin'};
            end
            
            simulationConditionIDs = repmat({simuConditionID}, [length(timepoints) 1]);
            
            obsPars_tmp = strsplit(ar.model(imodel).data(idata).fy{iy}, {'+','-','*','/','(',')','^',' '});
            obsPars_tmp = intersect(obsPars_tmp, ar.pLabel);
            obsPars_tmp = strcat(obsPars_tmp, ';');
            observableParameters = repmat({[obsPars_tmp{:}]}, [length(timepoints) 1]);
            
            noisePars_tmp = strsplit(ar.model(imodel).data(idata).fystd{iy}, {'+','-','*','/','(',')','^',' '});
            noisePars_tmp = intersect(noisePars_tmp, ar.pLabel);
            noisePars_tmp = strcat(noisePars_tmp, repmat(';', [length(noisePars_tmp)-1 1]));
            noiseParameters = repmat({[noisePars_tmp{:}]}, [length(timepoints) 1]);
            
            noiseDist = cell(length(timepoints),1);
            noiseDist(:) = {'normal'};
            
            measT = [measT; table(observableIDs,simulationConditionIDs,...
                measurements,timepoints,observableParameters, ...
                noiseParameters, obsTrafos, noiseDist)];
        end
    end
    
    condT = array2table(peCondValues);
    condT.Properties.VariableNames = peConds;
    
    condT = [table(conditionID'), condT];
    condT.Properties.VariableNames{1} = 'conditionId';
    
    measT.Properties.VariableNames = {'observableId', 'simulationConditionId', ...
        'measurement', 'time', 'observableParameters', 'noiseParameters',...
        'observableTransformation', 'noiseDistribution'};
    
    writetable(condT, ['SBML_Export/peTABcond_model' num2str(imodel) '.tsv'],...
        'Delimiter', '\t', 'FileType', 'text')
    writetable(measT, ['SBML_Export/peTABmeas_model' num2str(imodel) '.tsv'],...
        'Delimiter', '\t', 'FileType', 'text')
end
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

writetable(parT, 'SBML_Export/peTABpars_model.tsv',...
    'Delimiter', '\t', 'FileType', 'text')

%% Visualization Table
% for imodel = 1:length(ar.model)
%     for iplot = 1:length(ar.model(imodel).plot)
%         for idata = 1:length(ar.model(imodel).plot(iplot).dLink)
%             for iy = 1:length(ar.model(imodel).data(idata).y)
%                 plotID(end+1) =
%                 plotName(end+1) = [ar.model(imodel).data(idata).yNames{iy}, ...
%                     '_' ar.model(imodel).plot(iplot).condition{}];
%                 plotTypeSimulation(end+1) =
%                 plotTypeData(end+1) =
%                 datasetId(end+1) =
%                 independentVariable(end+1) =
%                 [independentVariableOffset] =
%                 [independentVariableName] =
%                 [legendEntry] = ar.model(imodel).plot(iplot).condition{};
%             end
%         end
%     end
% end
%
%writetable(visT, 'peTAB_visualization.tsv', 'Delimiter', '\t', 'FileType', 'text')

end