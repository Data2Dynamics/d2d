function out = measurement_table_has_timepoint_specific_mappings( ...
    measurement_df) %-> bool
    %Are there time-point or replicate specific assignments in the
    %measurement table.
    %
    %Arguments:
    %   measurement_df [table]:
    %       PEtab measurement table.
    %
    %Returns:
    %   bool
    %       True if there are time-point or replicate specific parameter
    %       assignments in the measurement table, false otherwise.
    
    cols = measurement_df.Properties.VariableNames;
    if ismember('noiseParameters', cols)    
        test = map(@isnumeric, measurement_df.noiseParameters);
        
        if iscell(measurement_df.noiseParameters)
            measurement_df.noiseParameters(test) = {NaN};
        elseif isnumeric(measurement_df.noiseParameters)
            measurement_df.noiseParameters(test) = NaN;
        end
    else
        measurement_df.noiseParameters = NaN(height(measurement_df), 1);
    end
    
    grouping_cols = get_notnull_columns(measurement_df, ...
        ["observableId", "simulationConditionId", ...
        "preequilibrationConditionId", "observableParameters", ...
        "noiseParameters"]);
    grouped_df = unique(measurement_df(:, grouping_cols), ...
        'rows', 'stable');
    
    grouping_cols = get_notnull_columns(grouped_df, ["observableId", ...
        "simulationConditionId", "preequilibrationConditionId"]);
    grouped_df2 = unique(grouped_df(:, grouping_cols), 'rows', ...
        'stable');
    
    if height(grouped_df) ~= height(grouped_df2)
        warning('Measurement table has timepoint-specific mappings');       
        out = true;
    else
        out = false;
    end
end