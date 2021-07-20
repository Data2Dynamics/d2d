function out = get_default_condition_file_name(model_name, folder) %-> string
    %Get file name according to proposed convention.
    %
    %Arguments:
    %   model_name string:
    %       Name of the model.
    %   folder string:
    %       Path of the folder that contains the model files.
    %
    %Returns:
    %   string:
    %       Path of the condition file.
    
    warning('GET_DEFAULT_CONDITION_FILE_NAME:DeprecationWarning', ...
        'This function will be removed in future releases')
    
    filename = sprintf('experimentalCondition_%s.tsv', model_name);
    out = fullfile(folder, filename);
end