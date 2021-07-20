function out = cvarsForCondition(problem, condid)
    %Returns condition dependent parameters for given condition id, as
    %defined in condition table.
    %
    %Arguments:
    %   problem Problem:
    %       PEtab problem object.
    %   condid string:
    %       Condition identifyer.
    %
    %Returns:
    %   [numeric]
    %       condition dependent parameter values.
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    x = string(problem.amimodel.x);
    
    cvars = setdiff(problem.condition_df.Properties.VariableNames, ...
        ["conditionId" "conditionName"]);
    cvars = setdiff(cvars, x, 'stable');
    
    condition_df = problem.condition_df;
    condition_df.Properties.RowNames = condition_df.conditionId;
    
    out = condition_df{condid, cvars};
end