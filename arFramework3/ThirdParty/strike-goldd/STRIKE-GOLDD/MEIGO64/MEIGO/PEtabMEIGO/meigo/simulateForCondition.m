function [t, x] = simulateForCondition(problem, condid, x0, p, is_preeq)
    %TODO
    
    if nargin == 4
        is_preeq = false;
    end
    
    opt_table = getOptimizationTableForCondition(problem, condid);
    cvars = cvarsForCondition(problem, condid);
    
    if is_preeq
        tode = [1E-8; 1E8];
    else
        texp = unique(opt_table.time, 'sorted');
        if numel(texp) == 1
            tode = [1E-8; texp];
        else
            tode = union(texp, 1E-8, 'sort');
        end
    end
    
    str = [problem.sbml_model.id '_dyn'];
    simfunc = @(t, x) feval(str, t, x, p, cvars);
    
    [t, x] = ode15s(simfunc, tode, x0);
    
    tdict = Dict();
    for i = 1:numel(t)
        tdict(string(t(i))) = x(i, :);
    end
    
    if is_preeq
        t = tode(end);
        x = x(end, :);
    else
        t = texp;
        
        n = size(x);
        x = zeros(numel(t), n(2));
        for i = 1:numel(t)
            x(i, :) = tdict(string(t(i)));
        end
    end
end