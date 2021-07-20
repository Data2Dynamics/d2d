function out = get_optimization_parameter_scales(parameter_df) %-> [Dict]
    %Get Dictionary with optimization parameter IDs mapped to parameter
    %scaling strings.
    %
    %Arguments:
    %   parameter_df [table]:
    %       PEtab parameter table.
    %
    %Returns:
    %   [Dict]
    %       Dictionary with optimization parameter IDs mapped to parameter
    %       scaling strings.
    
    parameter_scales = parameter_df.parameterScale( ...
        parameter_df.estimate == 1);
    
    out = Dict(get_optimization_parameters(parameter_df), ...
        parameter_scales);
end

