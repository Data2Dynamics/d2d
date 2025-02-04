function arSubstituteConditionsInModelFile(filename):
    %% Read file
    % Define all possible section titles
    SectionTitles = ["DESCRIPTION", "PREDICTOR", "COMPARTMENTS", "STATES", ...
                     "INPUTS", "REACTIONS", "DERIVED", "OBSERVABLES", ...
                     "ERRORS"];
    
    % Open the file
    fid = fopen(filename, 'r');
    
    % Check if file opened successfully
    if fid == -1
        error('Failed to open file');
    end
    
    % Initialize variables
    lines = {};
    linesConditions = {};
    conditionsSection = false;
    
    % Read file line by line
    while ~feof(fid)
        tline = fgetl(fid);
        lines{end+1} = tline;
    
        % Check if the line contains "CONDITION"
        if contains(tline, 'CONDITION')
            conditionsSection = true;
            continue; % Skip the CONDITION line itself
        end
    
        if contains(tline, SectionTitles)
            conditionsSection = false;
            continue;
        end
        
        % If we have found "CONDITION", start storing lines
        if conditionsSection
            parts = split(tline);
            if length(parts) >= 2
                if isnan(str2double(strrep(parts{2}, '"', ''))) % Check if number or equation is assigned to condition
                    linesConditions{end+1} = tline; % save equation for later
                    lines{end} = ''; % delete equation for new file
                end
            end
        end
    end
    
    % Close the file
    fclose(fid);
    
    %% Replace Conditions
    linesConditions = linesConditions(~strcmp(linesConditions, '        '));
    for i=1:length(linesConditions)
        condition = split(linesConditions{i});
        A = condition{1}; % name of contion
        B =  strrep(join(condition(2:end), ''), '"', ''); % Equation that will replace B
        B = ['(', B{1},')']; % add ( )
        lines = cellfun(@(x) strrep(x, A, B), lines, 'UniformOutput', false); % replace A with B in each line
    end
    
    %% write new file
    filename_def = split(filename,'.');
    filename_new = [filename_def{1},'_new.def'];
    
    % create
    fileID = fopen(filename_new, 'w'); 
    
    % Write each element of the cell array as a new line in the file
    for i = 1:length(lines)
        fprintf(fileID, '%s\n', lines{i});
    end
    
    % Close the file
    fclose(fileID);
