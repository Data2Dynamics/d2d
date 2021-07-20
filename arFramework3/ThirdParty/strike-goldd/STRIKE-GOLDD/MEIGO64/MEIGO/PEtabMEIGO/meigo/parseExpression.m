function out = parseExpression(problem, expression)
    %TODO
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    expression = [' ' char(expression) ' '];
    rex = '(?<=[\\W\\s])%s(?=[\\W\\s])';
    
    x = string(problem.amimodel.x);
    optpars = problem.x_free_ids;
    
    fixpars = problem.x_fixed_ids;
    fixvals = problem.x_nominal_fixed;
    
    comps = string({problem.sbml_model.compartment.id});
    compsize = ones(1, numel(comps));
    for i = 1:numel(comps)
        if problem.sbml_model.compartment(i).isSetSize
            compsize(i) = problem.sbml_model.compartment(i).size;
        end
    end
    
    cvars = setdiff(problem.condition_df.Properties.VariableNames, ...
        ["conditionId" "conditionName"], 'stable');
    cvars = setdiff(cvars, x, 'stable');
    
    [fixcomps, mask] = setdiff(comps, cvars, 'stable');
    fixcompsize = compsize(mask);
    
    for i = 1:numel(fixcompsize)
        expression = regexprep(expression, sprintf(rex, fixcomps(i)), ...
            sprintf('%f', fixcompsize(i)));
    end
    
    for i = 1:numel(x)
        expression = regexprep(expression, sprintf(rex, x(i)), ...
            sprintf('x(%d)', i));
    end
    
    for i = 1:numel(cvars)
        expression = regexprep(expression, sprintf(rex, cvars(i)), ...
            sprintf('c(%d)', i));
    end
    
    for i = 1:numel(fixpars)
        expression = regexprep(expression, sprintf(rex, fixpars(i)), ...
            sprintf('%f', fixvals(i)));
    end
    
    for i = 1:numel(optpars)
        expression = regexprep(expression, sprintf(rex, optpars(i)), ...
            sprintf('p(%d)', i));
    end    
    
    init_states_ids = x;
    for i = 1:numel(x)
        state = problem.sbml_model.getSpecies(x(i));        
        init_states_ids(i) = lower([char(x(i)) '0']);
        
        if state.isSetInitialConcentration
            tmp = state.initialConcentration;
        else
            tmp = 1;
        end
        
        expression = regexprep(expression, sprintf(rex, ...
            init_states_ids(i)), sprintf('%f', tmp));
    end
    
    otherpars = setdiff(problem.get_model_parameters, ...
        union(problem.x_ids, init_states_ids, 'stable'), 'stable');
    for i = 1:numel(otherpars)
        mask = strcmp({problem.sbml_model.parameter.id}, otherpars(i));
        
        if problem.sbml_model.parameter(mask).isSetValue
            tmp = problem.sbml_model.parameter(mask).value;
        else
            tmp = 1;
        end
        
        expression = regexprep(expression, sprintf(rex, ...
            otherpars(i)), sprintf('%f', tmp));
    end
    
    out = expression(2:end - 1);
end