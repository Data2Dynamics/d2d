function out = parsex0(problem, condid, preeq_x0)
    %TODO
    
    if nargin == 2
        preeq_x0 = [];
    end
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    if isempty_ext(preeq_x0)
        x0 = string(problem.amimodel.x0);
    else
        x0 = preeq_x0;
    end
    
    x = string(problem.amimodel.x);
    cvars = setdiff(problem.condition_df.Properties.VariableNames, ...
        ["conditionId" "conditionName"], 'stable');
    
    condition_df = problem.condition_df(:, cvars);
    condition_df.Properties.RowNames = problem.condition_df.conditionId;
    
    cvals = condition_df{condid, cvars};
    cdict = Dict(cvars, cvals);
    
    x0ids = intersect(x, cvars, 'stable');
    
    if ~isempty_ext(x0ids)
        for i = 1:numel(x0ids)
            x0(strcmp(x, x0ids(i))) = cdict(x0ids(i));
        end
    end
    
    for i = 1:numel(x0)
        x0(i) = parseExpression(problem, x0(i));
    end
    
    out = x0;
end