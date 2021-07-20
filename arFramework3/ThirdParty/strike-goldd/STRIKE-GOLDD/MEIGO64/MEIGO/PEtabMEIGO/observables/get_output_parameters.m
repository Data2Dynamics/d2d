function out = get_output_parameters(observable_df, sbml_model) %-> [string]
    %Get output parameters.
    %
    %Returns IDs of parameters used in observable and noise formulas not 
    %defined in the SBML model.
    %
    %Arguments:
    %   observable_df [table]:
    %       PEtab observable table.
    %   sbml_model [Sbml]:
    %       SBML model.
    %
    %Returns:
    %   [string]
    %       List of output parameter IDs.
    
    formulas = observable_df.observableFormula;
    if ismember('noiseFormula', observable_df.Properties.VariableNames)
        formulas = union(formulas, string(observable_df.noiseFormula));
    end
    
    tmp = unique(flatten(map(@symvar, formulas, 'UniformOutput', ...
        false)), 'stable');
    out = tmp(~map(@sbml_model.isElementById, tmp));    
end

