function out = get_parameter_df(parameter_file_name) %-> [table]
    %Read the provided parameter file into a 'table'
    %
    %Arguments:
    %   parameter_file_name string: 
    %       Name of the file to read from.
    %
    %Returns:
    %   [table]
    %       Parameter table
    
    out = parameter_file_name;
    if isempty_ext(parameter_file_name) || istable(parameter_file_name)
        return
    elseif ~exist(parameter_file_name, 'file')
        error('GET_PARAMETER_DF:FileNotFoundError', ...
            'No such file or directory')
    end
    
    out = readtable(out, ...
        'FileType', 'text', ...
        'ReadVariableNames', true, ...
        'PreserveVariableNames', true, ...
        'Delimiter', 'tab');
    
    columns = string(out.Properties.VariableNames);
    if ~ismember('parameterId', columns)
        error('GET_PARAMETER_DF:MandatoryFieldNotInTableError', ...
            'Parameters table missing mandatory field "parameterId"')
    end
    
    assert_no_leading_trailing_whitespace(columns, 'parameter')
end