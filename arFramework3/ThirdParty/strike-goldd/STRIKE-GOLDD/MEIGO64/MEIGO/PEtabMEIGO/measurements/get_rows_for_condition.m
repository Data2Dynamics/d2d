
function out = get_rows_for_condition(measurement_df, ...
    condition) %-> [table]
    %Extract rows in 'measurement_df' for 'condition' according
    %to 'preequilibrationConditionId' and 'simulationConditionId' in
    %'condition'.
    % 
    %Arguments:
    %   measurement_df [table]:
    %       PEtab measurement table.
    %   condition [table]:
    %       Table with single row and columns
    %       'preequilibrationConditionId' and 'simulationConditionId'.
    %       Or dictionary with those keys.
    % 
    %Returns:
    %   [table]
    %       The subselection of rows in 'measurement_df' for the condition
    %       'condition'.
    
    out = table();
    condition.Count = [];
    
    if istable(condition)
        condnames = string(condition.Properties.VariableNames);
        condid = string(condition{1, :});
    else
        condnames = condition.keys;
        condid = string(condition.values);
    end
    
    if ismember('simulationConditionId', condnames)
        simidx = strcmp(condnames, 'simulationConditionId');
        measurement_df = measurement_df(strcmp( ...
            measurement_df.simulationConditionId, condid{simidx}), :);
        
        out = measurement_df;
    end
    
    if ismember('preequilibrationConditionId', condnames)
        preeqidx = strcmp(condnames, 'preequilibrationConditionId');
        measurement_df = measurement_df(strcmp( ...
            measurement_df.preequilibrationConditionId, ...
            condid{preeqidx}), :);
        
        out = measurement_df;
    end
end