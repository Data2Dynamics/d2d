function out = get_condition_df(condition_file_name) %-> [table]
    %Read the provided condition file into a table.
    %
    %Conditions are rows, parameters are columns.
    %
    %Parameters:
    %   condition_file_name Union[string, table, []]:
    %                        File name of PEtab condition file.
    %
    %Returns:
    %   [table]: 
    %       Condition table.
    
    out = condition_file_name;
    if isempty_ext(condition_file_name) || istable(condition_file_name)       
        return
    elseif ~exist(condition_file_name, 'file')
        error('GET_CONDITION_DF:FileNotFoundError', ...
            'No such file')
    end        
    
    out = readtable(condition_file_name, ...
        'FileType', 'text', ...
        'ReadVariableNames', true, ...
        'PreserveVariableNames', true, ...
        'Delimiter', 'tab');
    
    columns = string(out.Properties.VariableNames);
    if ~ismember('conditionId', columns)
        error('GET_CONDITION_DF:MandatoryFieldNotInTableError', ...
            'Condition table missing mandatory field "conditionId"')
    end    
    
    assert_no_leading_trailing_whitespace(columns, 'condition') 
end