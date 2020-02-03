function arExportPEtab(name, export_SBML)
% arExportPEtab(name, export_SBML)
%
% Exports model, data, experimental condition and measurement details to 
% PEtab format (https://github.com/ICB-DCM/PEtab)
%
%   name          string that will be prepended to all exported files
%   export_SBML   export model SBML file as part of PEtab specification [true]

% Not implemented yet:
% - Visualization table
% - Parameter table: priorType, priorParameters: not yet specified!

global ar

if ~exist('export_SBML') || isempty(export_SBML)
    export_SBML = true;
end
if ~exist('name') || isempty(name)
    name = '';
end


%% Write Export Directory
if(~exist('./PEtab', 'dir'))
    mkdir('./PEtab')
end

%% Export SBML model
if export_SBML
    arExportSBML([],name);
end

%% Export data, conditions & parameters
for imodel = 1:length(ar.model)
    %% Condition and Measurement Table
    % Collect conditions
    condPar = {}; condVal = {}; conditionID = {};
    for idata = 1:length(ar.model(imodel).data)
        for icond = 1:length(ar.model(imodel).data(idata).condition)
            condPar{end+1} = ar.model(imodel).data(idata).condition(icond).parameter;
        end
%         % catch d2d fix for predictors that are not time
%         if ~(strcmp(ar.model(imodel).data(idata).t, 'time') || strcmp(ar.model(imodel).data(idata).t, 't'))
%             for it = 1:length(ar.model(imodel).data(idata).tExp)
%                 condPar{end+1} = [ar.model(imodel).data(idata).t...
%                     '_' num2str(it)];
%             end
%         end
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
            observableId = repmat(ar.model(imodel).data(idata).y(iy), ...
                [length(ar.model(imodel).data(idata).yExp(:,iy)) 1]);
            rowsToAdd = [table(observableId)];

            measurement = ar.model(imodel).data(idata).yExpRaw(:, iy);
            time = ar.model(imodel).data(idata).tExp;
            
            % skip if measurement contains only NaN
            if sum(isnan(measurement)) == numel(measurement)
                continue
            end

            % pre-equiblibration
            if isfield(ar, 'ss_conditions')
                preEquilibrationId = cell(length(time),1);
                preEquilibrationId(:) = ...
                    {['model' num2str(imodel) '_data' ...
                    num2str(ar.model(imodel).condition(ar.model.ss_condition.src).dLink)]};
                rowsToAdd = [rowsToAdd, table(preEquilibrationId)];
            end
            
            simulationConditionId = repmat({simuConditionID}, [length(time) 1]);
            rowsToAdd = [rowsToAdd, table(simulationConditionId)];

            rowsToAdd = [rowsToAdd, table(measurement)];
            rowsToAdd = [rowsToAdd, table(time)];
            
            % observable parameters
            obsPars_tmp = strsplit(ar.model(imodel).data(idata).fy{iy}, {'+','-','*','/','(',')','^',' '});
            obsPars_tmp = intersect(obsPars_tmp, ar.pLabel);
            obsPars_tmp = strcat(obsPars_tmp, ';');
            observableParameters = repmat({[obsPars_tmp{:}]}, [length(time) 1]);
            rowsToAdd = [rowsToAdd, table(observableParameters)];
            
            % noise parameters
            expErrors = ar.model(imodel).data(idata).yExpStdRaw(:,iy);
            if ar.config.fiterrors == -1
                if sum(isnan(expErrors)) > 0
                    error('arExportPEtab: Cannot use ar.config.fiterrors == -1 with NaN in exp errors')
                end
                noiseParameters = expErrors;
            else
                if ar.config.fiterrors == 0 && sum(isnan(expErrors)) == 0
                    noiseParameters = expErrors;
                else
                    noisePars_tmp = strsplit(ar.model(imodel).data(idata).fystd{iy}, {'+','-','*','/','(',')','^',' '});
                    noisePars_tmp = intersect(noisePars_tmp, ar.pLabel);
                    noisePars_tmp = strcat(noisePars_tmp, repmat(';', [length(noisePars_tmp)-1 1]));
                    noiseParameters = repmat({[noisePars_tmp{:}]}, [length(time) 1]);
                end
            end
            rowsToAdd = [rowsToAdd, table(noiseParameters)];

            % observable trafos
            observableTransformation = cell(length(time),1);
            observableTransformation(:) = {'log10'};
            if ~ar.model(imodel).data(idata).logfitting(iy)
                observableTransformation(:) = {'lin'};
            end
            rowsToAdd = [rowsToAdd, table(observableTransformation)];

            % noise dists
            noiseDistribution = cell(length(time),1);
            noiseDistribution(:) = {'normal'}; % others not possible in d2d
            rowsToAdd = [rowsToAdd, table(noiseDistribution)];

%             % experiment id
%             experimentId = cell(length(time),1);
%             experimentId(:) = {ar.model(imodel).data(idata).name};
%             rowsToAdd = [rowsToAdd, table(experimentId)];
% 
%             % idenpendent variable id
%             indVariableId = cell(length(time),1);
%             indVariableId(:) = {'time'};
%             rowsToAdd = [rowsToAdd, table(indVariableId)];
%             
%             % for dose response measurements:
%             if isfield(ar.model(imodel).data(idata), 'response_parameter') && ...
%                     ~isempty(ar.model(imodel).data(idata).response_parameter)
%                 indVariableId(:) = {ar.model(imodel).data(idata).response_parameter};
%             end
            measT = [measT; rowsToAdd];
        end
    end
    
    condT = array2table(peCondValues);
    condT.Properties.VariableNames = peConds;
    
    condT = [table(conditionID'), condT];
    condT.Properties.VariableNames{1} = 'conditionId';
        
    writetable(condT, ['PEtab/' name 'cond_model' num2str(imodel) '.tsv'],...
        'Delimiter', '\t', 'FileType', 'text')
    writetable(measT, ['PEtab/' name 'meas_model' num2str(imodel) '.tsv'],...
        'Delimiter', '\t', 'FileType', 'text')
end
%% Parameter Table
parameterScale_tmp = cell(1, length(ar.qLog10));
parameterScale_tmp(:) = {'lin'};
parameterScale_tmp(ar.qLog10 == 1) = {'log10'};

clear nominalValue_tmp
nominalValue_tmp = ar.p;
nominalValue_tmp = (1-ar.qLog10).*nominalValue_tmp + ar.qLog10.*10.^nominalValue_tmp;

lowerBound_tmp = ar.lb;
lowerBound_tmp = (1-ar.qLog10).*lowerBound_tmp + ar.qLog10.*10.^lowerBound_tmp;

upperBound_tmp = ar.ub;
upperBound_tmp = (1-ar.qLog10).*upperBound_tmp + ar.qLog10.*10.^upperBound_tmp;

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

writetable(parT, ['PEtab/' name 'pars_model.tsv'],...
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