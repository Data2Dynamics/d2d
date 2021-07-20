function out = get_optimization_parameters(parameter_df) %-> [string]
    %Get list of optimization parameter IDs from parameter table.
    %
    %Arguments:
    %   parameter_df [table]: 
    %       PEtab parameter table
    %
    %Returns:
    %   [string]
    %       Array of IDs of parameters selected for optimization.
    
    out = transpose(string(parameter_df.parameterId( ...
        parameter_df.estimate == 1)));
end
