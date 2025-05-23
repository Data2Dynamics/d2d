function dataFilenames = arWriteDataFromPEtab(yamlFile, savePath)
% arWriteDataFromPEtab(yamlFile, savePath)
%
% This function reads the PEtab files and writes the data to excel (.xls) files and creates the corresponding .def files for each simulation condition.
%
% Inputs:
%   yamlFile: string, name of the yaml file in the PEtab folder. If not specified, the function will look for a .yaml file in the PEtab folder. Multiple .yaml files are not allowed.
%   savePath: string, path to the folder where the data files should be saved. [default: 'DataPEtab']
%
% Outputs:
%   dataFilenames: cell array of strings, filenames of the data files
%  


%% Filenames
if ~exist('yamlFile','var') || isempty(yamlFile)
    yamlFile = dir(fullfile('PEtab', ['*.yaml']))
    if isempty(yamlFile)
        error('Did not find any .yaml file in PEtab folder.')
    elseif length(yamlFile)>1
        error('Did find multiple .yaml files in PEtab folder. Please specify!')
    end
else
    yamlFile = fullfile('PEtab', strrep(['PEtab',filesep,yamlFile],['PEtab',filesep],'')); % add PEtab folder if missing
    yamlFile = dir([strrep(yamlFile,'.yaml','') '.yaml']); % add .yaml file ending if missing
    if isempty(yamlFile)
        error('Did not find .yaml file.')
    end
end
yamlPath = yamlFile.folder;
yamlContent = arReadPEtabYaml(fullfile(yamlPath,yamlFile.name));

% check number of files per category, only one per category allowed atm
petabFiles = {'sbml_files', 'observable_files', 'measurement_files',...
    'condition_files', 'parameter_file'};
inputArgs = {};
for i = 1:numel(petabFiles)
    [out,numberOfEls] = extractFromStruct(yamlContent, petabFiles{i});
    if numberOfEls > 1
        error('arImportPEtab: YAML file can only refer to one file per category.')
    end
    inputArgs{end+1} = out;
end

measurementFile = yamlContent.measurement_files{1};
observableFile = yamlContent.observable_files{1};
conditionFile = yamlContent.condition_files{1};

% Path for d2d data files
if ~exist('savePath','var') || isempty(yamlFile)
    savePath = fullfile('DataPEtab');
end
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%% Model
m = TranslateSBML(fullfile(yamlPath,yamlContent.sbml_files{1}));
species = {};
for s=1:length(m.species)
    species{end+1} = m.species(s).name;
end

%% Data xls files

%% Read PETab files: measurements, conditions, observables 
Tdat = tdfread(fullfile(yamlPath,measurementFile)); % data / measurement file
fns = fieldnames(Tdat);
for i = 1:length(fns)
    if ischar(Tdat.(fns{i}))
        Tdat.(fns{i}) = regexprep(string(Tdat.(fns{i})),' ','');
    end
end

Tcond = tdfread(fullfile(yamlPath,conditionFile)); % condition File
fns = fieldnames(Tcond);
for i = 1:length(fns)
    if ischar(Tcond.(fns{i}))
        Tcond.(fns{i}) = regexprep(string(Tcond.(fns{i})), ' ', ''); 
    end
    if ismember(fns{i}, species)
            Tcond = renameStructField(Tcond, fns{i}, ['init_', fns{i}]);
    end
end

Tobs = tdfread(fullfile(yamlPath,observableFile)); % observable File
fns = fieldnames(Tobs);
for i = 1:length(fns)
    if ischar(Tobs.(fns{i}))
        Tobs.(fns{i}) = regexprep(string(Tobs.(fns{i})),' ','');
    end
end

%% For each unique simulationConditionId (second cloumn in table Tdat)
% create excel file:
uniqueSimConds = unique(Tdat.simulationConditionId); %Tcond.conditionId;
dataFilenames = {};
for i=1:length(uniqueSimConds)
    % if not existing create folder 'DataPEtab' in the current directory
    filename = fullfile(savePath, uniqueSimConds{i});
    dataFilenames{end+1} = filename;
    filenamexls = [filename, '.xls'];
    if exist(filenamexls, 'file')
        delete(filenamexls); % Delete the file
    end
    
    % find observabeles that of simulationConditionId
    observableId = unique(Tdat.observableId(contains(Tdat.simulationConditionId, uniqueSimConds{i})));

    ind = find(contains(Tdat.simulationConditionId, uniqueSimConds{i}));
    timeUni = unique(Tdat.time(ind));
    timeMax = [timeUni, NaN(length(timeUni),1)];
    
    for o = 1:length(observableId)
        indObs = intersect(find(Tdat.observableId == observableId{o}),ind); % observableID and simulationConditionId
        tObs = Tdat.time(indObs);
        
        % Compute unique values and counts
        [uniqueVals, ~, idx] = unique(tObs);
        counts = accumarray(idx, 1);
        
        % Fill counts into timeMax
        for k = 1:length(uniqueVals)
            matchIdx = timeMax(:,1) == uniqueVals(k); % Find matching rows in timeMax
            timeMax(matchIdx,2) = max(timeMax(matchIdx,2), counts(k)); % Assign count
        end
    end
    
    % time 
    time = repelem(timeMax(:,1), timeMax(:,2));
    [timeD1,timeD2] = size(time);
    if timeD2 ~= 1 
        time = time';
    end

    % create Table
    T = table(time, 'VariableNames', {'time'}); 
    for o = 1:length(observableId)
        T.(char(observableId(o))) = NaN(length(time), 1); % Fill with NaN for now
    end

    for o=1:length(observableId)
        indObs = intersect(find(Tdat.observableId == observableId{o}),ind); % observableID and simulationConditionId
        %indObs = find(Tdat.observableId == observableId{o});
        for t = 1:length(indObs)
        % find all rows in T where T.time==Tdat.time(indObs(t))
        indT = find(T.time == Tdat.time(indObs(t)));
            for tT = 1:length(indT)
                if isnan(T.(char(observableId(o)))(indT(tT)))
                    T.(char(observableId(o)))(indT(tT)) = Tdat.measurement(indObs(t));
                    break
                end
            end
        end
    end

    % writeTable to excel file
    writetable(T, filenamexls);

%% Data def files    
    if exist([filename, '.def'], 'file')
        delete([filename,'.def']); % Delete the file
    end
    fid = fopen([filename, '.def'], 'w');

    % Description
    fprintf(fid, '\nDESCRIPTION\n');
    fprintf(fid, '"This File is generated automatically from PEtab files"\n');

    % Predictor
    fprintf(fid, '\nPREDICTOR\n');
    %fprintf(fid, ['t\tT\tmin\ttime\t',num2str(min(time)),'\t',num2str(max(time)),'\n']);
    fprintf(fid, ['t\tT\tmin\ttime\t',num2str(0),'\t',num2str(max(time)*1.2),'\n']);

    % Inputs
    fprintf(fid, '\nINPUTS\n');

    % Observables   
    fprintf(fid, '\nOBSERVABLES\n');
    indObs = find(contains(Tobs.observableId, observableId));
    for j = 1:length(indObs)
        ind = indObs(j);
        comparelog10 = ~strcmp(Tobs.observableTransformation(ind),"lin");
        obsEntry = sprintf('%s\tC\tau\tconc.\t0\t    %u\t   \t"%s"', ...
            Tobs.observableId{ind},comparelog10,Tobs.observableFormula{ind});
        if isfield(Tobs,'observableName')
            fprintf(fid, [obsEntry,'\t"',Tobs.observableName{ind},'"\n']);
        else
            fprintf(fid, [obsEntry,'\n']);
        end
    end

    % Errors
    fprintf(fid, '\nERRORS\n');
    for j = 1:length(indObs)
        ind = indObs(j);
        if strcmp(Tobs.observableTransformation(ind),"log");
            fprintf(fid, [Tobs.observableId{ind},'\t','"',Tobs.noiseFormula{ind},' / log(10)','"\n']);
            warning(['Active fitting on natural logarithm scale for observable ',Tobs.observableId{ind},' in condition ',uniqueSimConds{i},'. Do not change ar.config.fiterrors!'])
        else
            if isnumeric(Tobs.noiseFormula)
                Tobs.noiseFormula = string(num2cell(Tobs.noiseFormula));
            end
            fprintf(fid, [Tobs.observableId{ind},'\t','"',Tobs.noiseFormula{ind},'"\n']);
        end
    end

    % Conditions            
    fprintf(fid, '\nCONDITIONS\n');
    indCond = find(contains(Tcond.conditionId, uniqueSimConds{i}));
    fns = fieldnames(Tcond);
    for fieldname = 2:length(fns)
        value = Tcond.(fns{fieldname})(indCond);
        if ((ischar(value) || isstring(value)) && strcmpi(value,"NaN"))
            value = str2num(value);
        end
        if ~(ischar(value) && isempty(value)) && ~(isnumeric(value) && isnan(value))
            fprintf(fid, [fns{fieldname}, '\t', '"',num2str(Tcond.(fns{fieldname})(indCond)),'"\n']);
        end
    end   
    fclose(fid);
end

end

function [out,numberOfEls] = extractFromStruct(struct, field)
if iscell(struct.(field))
    out = struct.(field);
    numberOfEls = numel(out);
else
    out = {struct.(field)};
    numberOfEls = 1;
end
end