function apply_output_parameter_overrides(mapping, cur_measurement_df)
    %Apply output parameter overrides to the parameter mapping dict for a 
    %given condition as defined in the measurement table 
    %("observableParameter","noiseParameters").
    %
    %Parameters:
    %   mapping (Dict): 
    %       parameter mapping dict as obtained from
    %       "get_parameter_mapping_for_condition".
    %   cur_measurement_df (table):
    %       Subset of the measurement table for the current condition.
        
    if ismember('observableParameters', ...
            cur_measurement_df.Properties.VariableNames)
        overrides = map(@split_parameter_replacement_list, ...
            cur_measurement_df.observableParameters, 'UniformOutput', false);
        wrapfunc = @(x, y) apply_overrides_for_observable(mapping, x, ...
            'observable', y);
        map(wrapfunc, cur_measurement_df.observableId, overrides);
    end
    
    if ismember('noiseParameters', ...
            cur_measurement_df.Properties.VariableNames)
        overrides = map(@split_parameter_replacement_list, ...
            cur_measurement_df.noiseParameters, 'UniformOutput', false);
        wrapfunc = @(x, y) apply_overrides_for_observable(mapping, x, ...
            'noise', y);
        map(wrapfunc, cur_measurement_df.observableId, overrides);
    end
end