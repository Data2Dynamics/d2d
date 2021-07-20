function out = get_default_parameter_file_name(model_name, folder) %-> string
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
    %       Path of the parameter file.
    
    warning('GET_DEFAULT_PARAMETER_FILE_NAME:DeprecationWarning', ...
        'This function will be removed in future releases')
    
    filename = sprintf('parameters_%s.tsv', model_name);
    out = fullfile(folder, filename);
end