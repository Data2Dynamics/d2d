function out = get_observable_df(observable_file_name) %-> [table]
    %Read the provided observable file into a MATLAB table.
    %
    %Arguments:
    %   observable_file_name string:
    %       Name of the file to read from.
    %
    %Returns:
    %   [table]
    %       Observable table.
    
    out = observable_file_name;
    if isempty_ext(observable_file_name) || istable(observable_file_name)
        return
    elseif ~exist(observable_file_name, 'file')
        error('GET_OBSERVABLE_DF:FileNotFoundError', ...
            'No such file or directory')
    end
    
    out = readtable(out, ...
        'FileType', 'text', ...
        'ReadVariableNames', true, ...
        'PreserveVariableNames', true, ...
        'Delimiter', 'tab');
    
    columns = string(out.Properties.VariableNames);
    if ~ismember('observableId', columns)
        error('GET_OBSERVABLE_DF:MandatoryFieldNotInTableError', ...
            'Observables table missing mandatory field "observableId"')
    end
    
    assert_no_leading_trailing_whitespace(columns, 'observable')
end