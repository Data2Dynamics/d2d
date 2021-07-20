function out = get_measurement_df(measurement_file_name) %-> [table]
    %Read the provided measurement file into a 'table'.
    %
    %Arguments:
    %   measurement_file_name [string, table, []]:
    %                          Name of file to read from.
    %
    %Returns:
    %   [table]
    %       Measurements table.
    
    out = measurement_file_name;
    if isempty_ext(measurement_file_name) || istable(measurement_file_name)
        return
    elseif ~exist(measurement_file_name, 'file')
        error('GET_MEASUREMENT_DF:FileNotFoundError', ...
            'No such file or directory')
    end
    
    out = readtable(measurement_file_name, ...
        'FileType', 'text', ...
        'ReadVariableNames', true, ...
        'PreserveVariableNames', true, ...
        'Delimiter', 'tab');
    
    columns = string(out.Properties.VariableNames);
    assert_no_leading_trailing_whitespace(columns, 'measurement')
end
