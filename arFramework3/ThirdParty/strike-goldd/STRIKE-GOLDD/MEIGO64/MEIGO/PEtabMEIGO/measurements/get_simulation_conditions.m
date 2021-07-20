function out = get_simulation_conditions(measurement_df) %-> [table]
    %Create a table of separate simulation conditions. A simulation
    %condition is a specific combination of simulationConditionId and
    %preequilibrationConditionId.
    %
    %Arguments:
    %   measurement_df [table]:
    %       PEtab measurement table.
    %
    %Returns:
    %   [table]
    %       Table with columns 'simulationConditionId' and
    %       'preequilibrationConditionId'. All null columns will be
    %       omitted.
    
    grouping_cols = get_notnull_columns(measurement_df, ...
        ["preequilibrationConditionId", "simulationConditionId"]);
    
    out = measurement_df(:, grouping_cols);
    [out, ~, ic] = unique(out, 'rows', 'stable');    
    
    out.Count = accumarray(ic, 1);
end