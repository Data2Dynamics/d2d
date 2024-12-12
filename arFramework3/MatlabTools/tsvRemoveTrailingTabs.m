function tsvRemoveTrailingTabs(inputFile, outputFile)
% tsvRemoveTrailingTabs removes trailing tabs (\t) from each line of a .tsv file
% and saves the cleaned data to a new file or overwrites the input file if
% no output file is specified.
%
% Usage:
%   tsvRemoveTrailingTabs('inputFile.tsv')               % Overwrites input file
%   tsvRemoveTrailingTabs('inputFile.tsv', 'outputFile.tsv') % Saves to output file
%
% Arguments:
%   inputFile  - The path to the input .tsv file.
%   outputFile - (Optional) The path to the output .tsv file. If not provided,
%                the input file is overwritten.

    % If outputFile is not specified, set it to inputFile (overwrite mode)
    if ~exist('outputFile','var') || isempty(outputFile)
        outputFile = inputFile;
    end

    % Open the input file for reading
    fid_in = fopen(inputFile, 'r');
    if fid_in == -1
        error('Could not open input file: %s', inputFile);
    end

    % Read all lines from the input file
    lines = {};
    while ~feof(fid_in)
        line = fgetl(fid_in);       % Read a line
        if ischar(line)            % Ensure the line is valid
            line = regexprep(line, '\t+$', ''); % Remove trailing tabs
            lines{end+1} = line;   % Add the cleaned line to the list
        end
    end
    fclose(fid_in); % Close the input file

    % Open the output file for writing
    fid_out = fopen(outputFile, 'w');
    if fid_out == -1
        error('Could not open output file: %s', outputFile);
    end

    % Write cleaned lines to the output file
    for i = 1:length(lines)
        fprintf(fid_out, '%s\n', lines{i});
    end
    fclose(fid_out);
end
