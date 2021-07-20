function apply_condition_parameters(par_mapping, scale_mapping, ...
    condition_id, condition_df, sbml_model)
    %Replace parameter IDs in parameter mapping dictionary by condition
    %table parameter values (in-place).
    %
    %Parameters:
    %   par_mapping (Dict): 
    %       see get_parameter_mapping_for_condition.
    %   condition_id (string, char): 
    %       ID of condition to work on.
    %   condition_df (table): 
    %       PEtab condition table.
    
    condition_df.Properties.RowNames = condition_df.conditionId;
    condition_df.conditionId = [];
    
    for i = 1:numel(condition_df.Properties.VariableNames)
        overridee_id = condition_df.Properties.VariableNames{i};
        
        if strcmp(overridee_id, 'conditionName')
            continue
        elseif ~isempty_ext(sbml_model.getSpecies(overridee_id))
            continue
        elseif ~isempty_ext(sbml_model.getCompartment(overridee_id))
            continue
        end
        
        par_mapping(overridee_id) = to_float_if_float(condition_df{ ...
            condition_id, overridee_id});
        scale_mapping(overridee_id) = 'lin';
    end
end