function y = arExportPEtab
% function y = arExportPEtab
%
% Exports data, condition and measurement details to PEtab format
% See also https://github.com/ICB-DCM/PEtab

%% Condition Table

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

% Write condition values
peCondValues = [];
numOfConds = length(peConds);
conditionID = {};

for imodel = 1:length(ar.model)
    for idata = 1:length(ar.model(imodel).data)
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
        observableNames = repmat(ar.model(imodel).data(idata).y, ...
            [length(ar.model(imodel).data(idata).yExp) 1]);
        observableNames = observableNames(:);
        simuConditionID
          
        ar.model(imodel).data(idata).yExp(:)
        ar.model(imodel).data(idata).tExp
        
        obsPars_tmp = cellfun(@strsplit, ar.model(imodel).data(idata).fy, repmat({{'+','-','*','/','(',')','^',' '}}, [length(ar.model(imodel).data(idata).fy) 1]), 'UniformOutput', false);
        obsPars_tmp = cellfun(@intersect, obsPars_tmp, repmat({ar.pLabel},[length(ar.model(imodel).data(idata).fy) 1]), repmat({'stable'},[length(ar.model(imodel).data(idata).fy) 1]), 'UniformOutput', false);
        obsPars = obsPars_tmp;
        
        noisePars_tmp = strsplit(ar.model(imodel).data(idata).fystd{iy}, {'+','-','*','/','(',')','^',' '});
        noisePars_tmp = intersect(noisePars_tmp, ar.pLabel);
        noisePars = noisePars_tmp;
        
        
        
    end
end

condT = array2table(peCondValues);
condT.Properties.VariableNames = peConds;

condT = [table(conditionID') condT];
condT.Properties.VariableNames{1} = 'conditionID';

writetable(condT, 'peTAB_conditions.tsv', 'Delimiter', '\t', 'FileType', 'text')


%% Condition Table 2
% 
% condT = table;
% 
% condPars = {}; condVals = {}; peConds = {};
% 
% % for idata = 1:length(ar.model.data)
% %     if ~isempty(ar.model.data(idata).condition)
% %         condPars{end+1} = ar.model.data(idata).condition.parameter;
% %         condVals{end+1} = ar.model.data(idata).condition.value;
% %     end
% % end
% 
% for imodel = 1:length(ar.model)
%     for idata = 1:length(ar.model(imodel).data)
%         if ~isempty(ar.model(imodel).data(idata).condition)
%             condPar = {ar.model(imodel).data(idata).condition.parameter};
%             condVal = {ar.model(imodel).data(idata).condition.value};
%             
%             condPars{end+1} = condPar;
%             condVals{end+1} = condVal;
%             
%             for icond = 1:length(condPar)
%                 if ~ismember(condPar{icond}, peConds)
%                     peConds{end+1} = condPar{icond};
%                 end
%                 condPos = find(contains(peConds, condPar{icond}));
%             end
%             
%             
%             
%         end
%     end
% end
% 
% 
% for imodel = 1:length(ar.model)
%     for idata = 1:length(ar.model(imodel).data)
%         if ~isempty(ar.model(imodel).data(idata).condition)
%             
%         end
%     end
% end

%% Table of Measurements
for imodel = 1:length(ar.model)
    for idata = 1:length(ar.model(imodel).data)
        for icond = 1:length(ar.model(imodel).data(idata).condition)
            for iy = 1:length(ar.model(imodel).data(idata).y)
                
                observableID = ar.model(imodel).data(idata).y{iy};
                simulationConditionId = [ar.model(imodel).data(idata).condition(icond).parameter ...
                    '_' ar.model(imodel).data(idata).condition(icond).value];
                measurement = ar.model(imodel).data(idata).yExp(:,iy);
                time = ar.model(imodel).data(idata).tExp;
                
                obsPars_tmp = strsplit(ar.model(imodel).data(idata).fy{iy}, {'+','-','*','/','(',')','^',' '});
                obsPars_tmp = intersect(obsPars_tmp, ar.pLabel);
                obsPars = obsPars_tmp;

                noisePars_tmp = strsplit(ar.model(imodel).data(idata).fystd{iy}, {'+','-','*','/','(',')','^',' '});
                noisePars_tmp = intersect(noisePars_tmp, ar.pLabel);
                noisePars = noisePars_tmp;
                  
            end
        end
    end
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
    'parameterScale', 'lowerBound', 'upperBound', 'nominalValue', 'estimate'};

writetable(parT, 'peTAB_parameters.tsv', 'Delimiter', '\t', 'FileType', 'text')

end