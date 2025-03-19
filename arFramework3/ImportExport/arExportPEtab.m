function arExportPEtab(name, export_SBML)
% arExportPEtab(name, export_SBML)
% Export parameter estimation problem to PEtab standard.
%
%   name          string that will be prepended to all exported files
%   export_SBML   export model SBML file [true]
%
% See also
%       arImportPEtab
%
% References
%   - PEtab standard: https://petab.readthedocs.io/en/latest/
%   - D2D Wiki: https://github.com/Data2Dynamics/d2d/wiki/Support-for-PEtab

global ar

if ~exist('export_SBML','var') || isempty(export_SBML)
    export_SBML = true;
end

if ~exist('name','var') || isempty(name)
    directory = split(pwd,filesep);
    name = directory{end};
end

%% Write Export Directory
if(~exist('./PEtab', 'dir'))
    mkdir('./PEtab')
end

%% Export SBML model
if export_SBML
    arExportSBML(name);
end

%% Models naming convention
IDs = struct();
IDs.name = name;
if length(ar.model) > 1
    IDs.model = cellfun(@(x) [name '_' x], {ar.model.name}, 'UniformOutput', false);
else
    IDs.model = {name};
end

%% Export observables, conditions, measurements and parameter table
for m = 1:length(ar.model)
    IDs = writeConditionsTable(m, IDs);
    IDs = writeObservablesTable(m, IDs);
    IDs = writeMeasurementsTable(m, IDs);
end
writeParameterTable(IDs);

%% Write YAML file
writeYAMLfile(IDs);

end


function IDs = writeConditionsTable(m, IDs)

global ar

% generate condition IDs:
% combine datafile name and condition number (needed for unique mapping)
numDigits = max(1, floor(log10(length(ar.model(m).data))) + 1);  % Ensure at least 1 digit
IDs.condition = arrayfun(@(d) sprintf(['D%0', num2str(numDigits), 'd_%s'], d, ar.model(m).data(d).name), ...
    1:length(ar.model(m).data), "UniformOutput", false);

% collect all data-condition-specific parameter replacements
condT = table(ar.model(m).p, ar.model(m).fp', 'VariableNames', {'modelP', 'modelFP'});
for d = 1:length(ar.model(m).data)
    is_modelP = ismember(ar.model(m).data(d).pold, ar.model(m).p);
    condT_tmp = cell2table([ar.model(m).data(d).pold(is_modelP)', ar.model(m).data(d).fp(is_modelP)]);
    condT_tmp.Properties.VariableNames = {'modelP', IDs.condition{d}};
    condT = outerjoin(condT, condT_tmp, "MergeKeys", true);
end

% replace init parameters with the corresponding stateID
for jp = 1:size(condT, 1)
    findInit = strcmp(condT.modelP(jp), ar.model(m).px0);
    if any(findInit)
        condT.modelP(jp) = ar.model(m).x(findInit);
    end
end

% remove lines where all replacements are identical to the global replacement
qRemove = arrayfun(@(jp) all(strcmp(condT.modelFP(jp), condT{jp, 3:end})), 1:size(condT, 1));
condT = condT(~qRemove, :);

% evaluate expressions (sometimes expressions reduce) and try to remove again
condT = array2table(arrayfun(@(x) string(arSym(x)), condT{:,:}), "VariableNames", condT.Properties.VariableNames);
qRemove = arrayfun(@(jp) all(strcmp(condT.modelFP(jp), condT{jp, 3:end})), 1:size(condT, 1));
condT = condT(~qRemove, :);

% reformat to match PEtab requirements
finalCondTab = array2table(transpose(condT{:,:}));
finalCondTab = [condT.Properties.VariableNames', finalCondTab];
if all(qRemove)
    finalCondTab.Properties.VariableNames = {'conditionId'};
else
    finalCondTab.Properties.VariableNames = ['conditionId' finalCondTab{1,2:end}];
end
finalCondTab = finalCondTab(3:end, :);

% export conditions table
filename = ['PEtab' filesep  IDs.model{m}  '_conditions.tsv'];
writetable(finalCondTab, filename, 'Delimiter', '\t', 'FileType', 'text')
disp([filename, ' written for model ', num2str(m)])

end


function IDs = writeObservablesTable(m, IDs)

global ar

IDs.obs = {};

obsT_tmp = table;

for d = 1:length(ar.model(m).data)
    % generate observable IDs:
    obsId = ar.model(m).data(d).yNames';
    IDs.obs{d} = cell(size(ar.model(m).data(d).yNames')); % for later use
    obsName = ar.model(m).data(d).yNames';

    obsFormula = ar.model(m).data(d).fy;
    for ify = 1:size(obsFormula,1)
        % replace parameter substitutions
        obsFormula{ify} = char(arSubs(arSym(obsFormula{ify}), ...
            arSym(ar.model(m).data(d).pold), ...
            arSym(ar.model(m).data(d).fp')));
        % replace derived quantities
        if ~isempty(ar.model(m).z)
            symFz = arSym(obsFormula{ify});
            symFzSubs = arSubs(symFz, arSym(ar.model(m).z), arSym(ar.model(m).fz'));
            obsFormula{ify} = char(symFzSubs);
        end
    end
    if size(obsFormula,2) > 1; obsFormula = obsFormula'; end

    obsTrafo = cell(length(obsId),1);
    obsTrafo(:) = {'lin'};
    obsTrafo(logical(ar.model(m).data(d).logfitting)) = {'log10'};

    noiseFormula = ar.model(m).data(d).fystd;
    for ify = 1:size(obsFormula,1)
        % replace parameter substitutions
        noiseFormula{ify} = char(arSubs(arSym(noiseFormula{ify}), ...
            arSym(ar.model(m).data(d).pold), ...
            arSym(ar.model(m).data(d).fp')));
        % replace derived quantities
        if ~isempty(ar.model(m).z)
            symFz = arSym(noiseFormula{ify});
            symFzSubs = arSubs(symFz, arSym(ar.model(m).z), arSym(ar.model(m).fz'));
            noiseFormula{ify} = char(symFzSubs);
        end
    end
    noiseFormulaSubs = cell(size(noiseFormula));
    for ifystd = 1:length(noiseFormula)
        noiseFormulaSymSingle = arSym(noiseFormula{ifystd});
        noiseFormulaSubs{ifystd} = char(arSubs(noiseFormulaSymSingle, arSym(obsName), arSym(obsId)));
    end
    noiseFormula = noiseFormulaSubs;
    if size(noiseFormula,2) > 1; noiseFormula = noiseFormula'; end

    noiseDistribution = cell(length(obsId),1);
    noiseDistribution(:) = {'normal'}; % others not possible in d2d

    obsT_tmp = [obsT_tmp; ...
        table(obsId, obsName, obsFormula, obsTrafo, ...
        noiseFormula, noiseDistribution)];
end

[obsT, ind_tmpA, ind_tmpC] = unique(obsT_tmp, 'rows');

% enumerate obsIds that occur more than once
[~, indA, indC] = unique(obsT.obsId, 'stable');
counts = accumarray(indC, 1);
occurrences = zeros(size(obsT.obsId));
for i = 1:height(obsT)
    if counts(indC(i)) > 1
        occurrences(i) = sum(indC(1:i) == indC(i));
        numDigits = floor(log10(counts(indC(i)))) + 1;
        obsT.obsId{i} = sprintf('%s_%0*d', obsT.obsId{i}, numDigits, occurrences(i));
    end
end

% Filling IDs.obs for later use in measurements
IDsObs_renamed = obsT.obsId(ind_tmpC, :);
count = 1;
for d=1:length(IDs.obs)
    for i=1:length(IDs.obs{d})
        IDs.obs{d}{i,:} = IDsObs_renamed{count};
        count = count+1;
    end
end

% write table to file
obsT.Properties.VariableNames = { ...
    'observableId', 'observableName', ...
    'observableFormula', 'observableTransformation', ...
    'noiseFormula', 'noiseDistribution'};
filename = ['PEtab' filesep IDs.model{m} '_observables.tsv'];
writetable(obsT, filename, 'Delimiter', '\t', 'FileType', 'text')
disp([filename, ' written for model ', num2str(m)])
end


function IDs = writeMeasurementsTable(m, IDs)

global ar

persistent threwNormWarning
if isempty(threwNormWarning)
    threwNormWarning = 0;
end

%% Measurement Table
measT = table;

for d = 1:length(ar.model(m).data)
    for iy = 1:length(ar.model(m).data(d).y)

        % observable ID
        time = ar.model(m).data(d).tExp;
        rowsToAdd = [table(repmat(IDs.obs{d}(iy), [length(time),1]), 'VariableNames', {'observableId'})];

        % convert measurements if necessary (log10, normalize)
        if ar.model(m).data(d).logplotting(iy)
            measurement = 10.^(ar.model(m).data(d).yExp(:, iy));
        else
            measurement = ar.model(m).data(d).yExp(:, iy);
        end
        % ar.model.yExp is already normalized, if normalization is active. Doing it again causes problems.
        % if ar.model(m).data(d).normalize(iy) == 1
        %     if ~threwNormWarning
        %         warning('Normalization of experimental measurements is not supported in PEtab. Measurement values in ar.model(:).data(:).yExpRaw will be normalized before export.')
        %         threwNormWarning = 1;
        %     end
        %     measurement = measurement/max(measurement);
        % end

        % skip if measurement contains only NaN
        if sum(isnan(measurement)) == numel(measurement)
            continue
        end

        % pre-equiblibration
        if isfield(ar, 'ss_conditions') && ar.ss_conditions
            for iss=1:length(ar.model(m).ss_condition)
                ss_d = ar.model(m).condition(ar.model(m).ss_condition(iss).src).dLink;
                preequilibrationConditionId = repmat(IDs.condition(ss_d), [length(time) 1]);
                rowsToAdd = [rowsToAdd, table(preequilibrationConditionId)];
            end
        end

        simulationConditionId = repmat(IDs.condition(d), [length(time) 1]);
        rowsToAdd = [rowsToAdd, table(simulationConditionId)];
        rowsToAdd = [rowsToAdd, table(measurement)];
        rowsToAdd = [rowsToAdd, table(time)];

        % observable parameters
        observableParameters = cell([length(time) 1]);
        rowsToAdd = [rowsToAdd, table(observableParameters)];

        % noise parameters
        expErrors = ar.model(m).data(d).yExpStd(:,iy);
        if ar.config.fiterrors == -1
            if sum(isnan(expErrors(~isnan(measurement)))) > 0
                error('arExportPEtab: Cannot use ar.config.fiterrors == -1 with NaN in exp errors')
            end
            noiseParameters = expErrors;
        else
            if ar.config.fiterrors == 0 && sum(isnan(expErrors)) ~= length(expErrors)
                noiseParameters = expErrors;
            else
                noiseParameters = cell([length(time) 1]);
            end
        end
        block = table(num2cell(noiseParameters));
        block.Properties.VariableNames = {'noiseParameters'};
        rowsToAdd = [rowsToAdd, block];

        % noise distributions
        noiseDistribution = cell(length(time),1);
        noiseDistribution(:) = {'normal'}; % others not possible in d2d
        rowsToAdd = [rowsToAdd, table(noiseDistribution)];

        % add to measurement table
        measT = [measT; rowsToAdd];

    end
end
filename = ['PEtab' filesep  IDs.model{m}  '_measurements.tsv'];
writetable(measT, filename, 'Delimiter', '\t', 'FileType', 'text')
disp([filename, ' written for model ', num2str(m)])

end


function writeParameterTable(IDs)
%% Parameter Table (one for the whole project, i.e. for all models jointly)

global ar

% paramter ID and name
parameterId = ar.pLabel;
parameterName = ar.pLabel;

% parameter flag for log10 trafo
parameterScale = cell(1, length(ar.p));
parameterScale(:) = {'lin'};
parameterScale(ar.qLog10 == 1) = {'log10'};

% parameter flag for fitting
estimate = ar.qFit;
estimate(ar.qFit == 2) = 0;

% value, lower and upper bounds are on linear scale
nominalValue = ar.p;
lowerBound = ar.lb;
upperBound = ar.ub;
nominalValue(ar.qLog10 == 1) = 10.^ar.p(ar.qLog10 == 1);
lowerBound(ar.qLog10 == 1) = 10.^ar.lb(ar.qLog10 == 1);
upperBound(ar.qLog10 == 1) = 10.^ar.ub(ar.qLog10 == 1);

% priors
% map d2d prior definations 0-3 to corresponding PEtab prior definations
priorTypes = {'parameterScaleUniform', 'normal', 'uniform', 'laplace'};
initializationPriorType = priorTypes(ar.type + 1);
% initializationPriorParameters
initializationPriorParameters = cell(size(ar.type));
for i = 1:length(ar.type)
    switch ar.type(i)
        case 0  % No prior
            initializationPriorParameters{i} = '';
        case 1  % Normal
            initializationPriorParameters{i} = sprintf('%f, %f', ar.mean(i), ar.std(i));
        case 2  % Uniform
            initializationPriorParameters{i} = sprintf('%f, %f', ar.lb(i), ar.ub(i));
        case 3  % Laplace
            initializationPriorParameters{i} = sprintf('%f, %f', ar.mean(i), ar.std(i));
    end
end

% build table and write to disk
variableNames = {'parameterId', 'parameterName', 'parameterScale', ...
    'lowerBound', 'upperBound', 'nominalValue', 'estimate', ...
    'initializationPriorType', 'initializationPriorParameters'};
parT = table(parameterId(:), parameterName(:), parameterScale(:), ...
    lowerBound(:), upperBound(:), nominalValue(:), estimate(:), ...
    initializationPriorType(:), initializationPriorParameters(:),...
    'VariableNames', variableNames);
filename = ['PEtab' filesep IDs.name '_parameters.tsv'];
writetable(parT, filename, 'Delimiter', '\t', 'FileType', 'text')
disp([filename, ' written'])

end

% function writeVisualisationTable(m, IDs)
% %% Visualisation Table
% % this is a very simple implementation, only covers the basic features like:
% % - combining datasets
% % - lin/log10 trafo
% % - time-course / dose-response

% % plotId, plotName, datasetId, xValues, xLabel, xScale, yValues, yLabel, yScale

% global ar

% for jp = 1:length(ar.model(m).plot)
%     plotId = sprintf('M%i_Plot%i', m, jp);
%     plotName = ar.model(m).plot(jp).name;

%     dLink = ar.model(m).plot(jp).dLink;
%     for jd = 1:length(dLink)
%         d = dLink(jd);
%         s
%     end


% end



% end


function writeYAMLfile(IDs)
%% YAML File
% one parameter table in total
% each "ar.model" is equivalent to one PEtab "problem"
% one SBML file for each PEtab problem
% one obs table for each PEtab problem
% one cond table (each condition coresp. ar.model.data) for each PEtab problem
% one meas table for each PEtab problem

global ar

% content of yaml file (line by line)
yamlLines = {};
yamlLines{end+1} = ['format_version: 1'];
yamlLines{end+1} = ['parameter_file: ' IDs.name '_parameters.tsv'];
yamlLines{end+1} = ['problems:'];
for m = 1:length(ar.model)
    yamlLines{end+1} = ['  - sbml_files:'];
    yamlLines{end+1} = ['    - ' IDs.model{m} '_model.xml'];
    yamlLines{end+1} = ['    observable_files:'];
    yamlLines{end+1} = ['    - ' IDs.model{m} '_observables.tsv'];
    yamlLines{end+1} = ['    condition_files:'];
    yamlLines{end+1} = ['    - ' IDs.model{m} '_conditions.tsv'];
    yamlLines{end+1} = ['    measurement_files:'];
    yamlLines{end+1} = ['    - ' IDs.model{m} '_measurements.tsv'];
end

% write yaml file to disk
filename = ['PEtab' filesep IDs.name '.yaml'];
fid = fopen(filename, 'w');
fprintf(fid, strjoin(yamlLines, '\n'));
fclose(fid);
disp([filename, ' written'])

end