function arRenameModelCondPars(m, targetDir, makeCopies)
% ARRENAMEMODELCONDITIONPARS removes dependencies between condition replacements
% by renaming the model parameters in the model.def (and if necessary data.def) file(s).
%
% Explanation:
% In d2d the replacements specified in CONDITIONS are applied exactly once and independently
% of each other. Let's illustrate this with the example of swapping two parameters in the model:
%   myExpression = str2sym('x^y');
%   p = arSym({'x', 'y'});
%   fp = arSym({'y', 'x'});
%   arSubs(myExpression, p, fp) returns 'y^x'
%
% However, if we want to export the model as a SBML file, the conditions can only be encoded as
% "initialAssignment" rules. These rules are all combined to form a set of equations that are
% solved simultaneously. Therefore, the previous example cannot be encoded in SBML since there
% is a circular dependency: x depends on y and y depends on x.
% 
% To resolve this, we want to rename the original parameters in the reactions, ODEs and inputs.
% In the example this would mean:
%   myExpression = str2sym('x_model^y_model');
%   p = arSym({'x_model', 'y_model'});
%   fp = arSym({'y', 'x'});
%   arSubs(myExpression, p, fp) returns 'y^x'
%
% Finally, since there might also be conditions on the model parameters in the data.def files,
% we also need to rename the parameters in the data.def files.

arguments
    m (1,1) double
    targetDir (1,1) string = ''
    makeCopies (1,1) logical = true
end


global ar


% Check if model is SBML compatible
[qSBMLCompatible, pNeedsRenaming] = arCheckSBMLCompatibility(m);
if qSBMLCompatible
    fprintf('Model conditions already SBML compatible. No parameters need renaming.\n')
    return
end

%% Overwrite the project or create a new one
if isempty(targetDir)
    % overwrite the project
    targetDir = ar.info.path;
else
    % copy the project to a new directory
    targetDir = fullfile(targetDir);
    if exist(targetDir, 'dir')
        warning('Directory already exists: %s\n', targetDir)
    end
    copyfile(ar.info.path, targetDir)
    makeCopies = false;
end


%% Read old and write new *.def files

% find correct model.def file
callLoadModel = strcmp(ar.setup.commands, 'arLoadModel');
if sum(callLoadModel) == 1
    modelDef = fullfile(ar.info.path, ar.setup.modelfiles{callLoadModel});
elseif sum(callLoadModel) > 1
    % multiple models loades, search for the correct one
    error('Multiple models loaded. Seraching the correct one is not implemented yet.')
else
    error('No model.def file found in ar.setup.modelfiles.')
end

% find all data.def files
callLoadData = strcmp(ar.setup.commands, 'arLoadData');
dataFiles = ar.setup.datafiles(callLoadData);
dataDefs = cellfun(@(d) fullfile(ar.info.path, d{1}), dataFiles, 'UniformOutput', false);
dataDefs = unique(dataDefs);  % remove duplicates (if multiple model.data share a data.def)
defFiles = [modelDef, dataDefs];

% flag for each parameter (did it appear in reactions, ODEs or inputs?)
qInModelReactions = false(1, length(pNeedsRenaming));
qInDataReactions = false(length(dataDefs), length(pNeedsRenaming));

for fileIdx = 1:length(defFiles)

    didReplace = false;
    defFile = defFiles{fileIdx};

    % backup and target file
    defFileTarget = strrep(defFile, ar.info.path, targetDir);
    defFileCopy = sprintf('%s_original.def', defFile(1:end-4));
    defFileCopy = strrep(defFileCopy, ar.info.path, targetDir);
    if makeCopies
        copyfile(defFile, defFileCopy);
    end

    % read file and split content into sections
    fileContent = fileread(defFile);
    defFileHeads = ["CONDITIONS", "PARAMETERS"];
    [sections, headings] = split(fileContent, defFileHeads);

    % modify file content 
    for ip = 1:length(pNeedsRenaming)
        paramOld = pNeedsRenaming{ip};
        paramNew = [paramOld, '_model'];

        % modify the first section (before CONDITIONS)
        % -> replace all occurrences of p by p_model
        pattern1 = sprintf('\\<%s\\>', paramOld);
        matches = regexp(sections{1}, pattern1, 'match');
        if ~isempty(matches)
            didReplace = true;
            if fileIdx == 1
                qInModelReactions(ip) = true;
            else
                qInDataReactions(fileIdx-1, ip) = true;
            end
            sections{1} = regexprep(sections{1}, pattern1, paramNew);
        end

        % modify CONDITIONS section (line by line backwards)
        % -> replace all occurrences of p by p_model in LHS of conditions
        % -> also keep old conditions if p is a state init parameter
        condSection = splitlines(sections{2});
        pattern2 = sprintf('^\\<%s\\>', paramOld);
        for il = length(condSection):-1:1
            oldCond = condSection{il};
            matches = regexp(oldCond, pattern2, 'match');
            if ~isempty(matches)
                didReplace = true;
                if ~ismember(paramOld, ar.model(m).px0)
                    % normal parameter
                    % -> replace p by p_model in LHS of conditions
                    condSection{il} = regexprep(oldCond, pattern2, paramNew);
                elseif qInModelReactions(ip) || qInDataReactions(fileIdx-1, ip)
                    % init param that also appears in reactions
                    % -> keep the old condition (for initial value of state)
                    % -> use p_model as new parameter in equations
                    condSection{il} = regexprep(oldCond, pattern2, paramNew);
                    condSection = [condSection(1:il); oldCond; condSection(il+1:end)];
                else
                    % init param that does not appear in reactions
                    % -> do nothing, no renaming necessary
                end
            end
        end
        sections(2) = join(condSection, newline);
    end

    % console output
    [~, fileName] = fileparts(defFile);
    if didReplace
        fprintf('%s.def: did replace model parameters.\n', fileName);
    else
        fprintf('%s.def: no replacements. No backup necessary.\n', fileName);
        if makeCopies
            delete(defFileCopy)
        end
    end

    % re-join the sections and write to the original file#
    fileContent = join(sections, headings);
    fid = fopen(defFileTarget, 'w');
    fprintf(fid, '%s', fileContent{:});
    fclose(fid);

end


%% Read and modify data tables (if necessary)
dataTables = cellfun(@(d) fullfile(ar.info.path, d{2}), dataFiles, 'UniformOutput', false);
dataTables = unique(dataTables);  % remove duplicates (if multiple model.data share a data table)

for d = 1:length(dataTables)
    
    dataFile = dataTables{d};
    [~, dataName, ext] = fileparts(dataFile);

    % backup and target file
    dataFileCopy = sprintf('%s_original%s', dataFile(1:end-4), ext);
    dataFileCopy = strrep(dataFileCopy, ar.info.path, targetDir);
    dataFileTarget = strrep(dataFile, ar.info.path, targetDir);

    % read, modify and write data file
    didReplace = false;
    switch ext
        case {'.xls', '.xlsx'}
            [~, ~, data] = xlsread(dataFile);
            for ip = 1:length(pNeedsRenaming)
                paramOld = pNeedsRenaming{ip};
                paramNew = [paramOld, '_model'];
                qParam = strcmp(data(1,:), paramOld);
                if any(qParam)
                    didReplace = true;
                    if ~ismember(paramOld, ar.model(m).px0)
                        % normal parameter
                        % -> replace p by p_model in column names
                        data{1, qParam} = paramNew;
                    elseif qInModelReactions(ip) || qInDataReactions(fileIdx-1, ip)
                        % init param that also appears in reactions
                        % -> keep the old column (for initial value of state)
                        % -> use p_model as new parameter in equations
                        idReplace = find(qParam);
                        dataInsert = data(:, idReplace);
                        dataInsert{1} = paramNew;
                        data = [data(:, 1:idReplace), dataInsert, data(:, idReplace+1:end)];
                    else
                        % init param that does not appear in reactions
                        % -> do nothing, no renaming necessary
                    end
                end
            end
            if didReplace
                if makeCopies
                    copyfile(dataFile, dataFileCopy);
                end
                xlswrite(dataFileTarget, data);
            end
        case '.csv'
            data = readtable(dataFile);
            for ip = 1:length(pNeedsRenaming)
                paramOld = pNeedsRenaming{ip};
                paramNew = [paramOld, '_model'];
                varNames = data.Properties.VariableNames;
                qParam = strcmp(varNames, paramOld);
                if any(qParam)
                    didReplace = true;
                    if ~ismember(paramOld, ar.model(m).px0)
                        % normal parameter
                        % -> rename column
                        varNames{qParam} = paramNew;
                        data.Properties.VariableNames = varNames;
                    elseif qInModelReactions(ip) || qInDataReactions(fileIdx-1, ip)
                        % init param that also appears in reactions
                        % -> keep the old column (for initial value of state)
                        % -> use p_model as new parameter in equations
                        % -> add duplicate column with new name
                        idReplace = find(qParam);
                        dataInsert = data(:, idReplace);
                        dataInsert.Properties.VariableNames = {paramNew};
                        data = [data(:, 1:idReplace), dataInsert, data(:, idReplace+1:end)];
                    else
                        % init param that does not appear in reactions
                        % -> do nothing, no renaming necessary
                    end
                end
            end
            if didReplace
                if makeCopies
                    copyfile(dataFile, dataFileCopy);
                end
                writetable(data, dataFile);
            end
    end
    if didReplace
        fprintf('%s%s: did replace model parameters in column names.\n', dataName, ext);
    else
        fprintf('%s%s: no replacements. No backup necessary.\n', dataName, ext);
    end

end
