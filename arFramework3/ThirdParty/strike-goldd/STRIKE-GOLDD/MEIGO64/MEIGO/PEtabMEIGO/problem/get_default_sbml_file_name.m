function out = get_default_sbml_file_name(model_name, folder) %-> string
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
    %       Path of the sbml model.
    
    warning('GET_DEFAULT_SBML_FILE_NAME:DeprecationWarning', ...
        'This function will be removed in future releases')
    
    filename = sprintf('model_%s.xml', model_name);
    out = fullfile(folder, filename);
end