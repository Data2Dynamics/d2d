function out = parseObservableFormulas(problem, opt_problem, condid)
    %TODO
    
    out = opt_problem;
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    rex = '(?<=[\\W\\s])%s(?=[\\W\\s])';
    
    condition_df = problem.condition_df;
    condition_df.Properties.RowNames = condition_df.conditionId;
    condition_df.conditionId = [];
    
    cvars = setdiff(problem.condition_df.Properties.VariableNames, ...
        ["conditionId", "conditionName"], 'stable');
    cvars = Dict(cvars, condition_df{condid, cvars});
    
    for i = 1:height(out)
        obsid = out.observableId{i};
        formula = [' ' out.observableFormula{i} ' '];
        
        if ismember('observableParameters', out.Properties.VariableNames)
            if isnumeric(out.observableParameters)
                obspar = out.observableParameters(i);
            else
                obspar = out.observableParameters{i};
            end
            obspar = split_parameter_replacement_list(obspar);          
            
            for j = 1:numel(obspar)
                overridee = sprintf('observableParameter%d_%s', j, obsid);
                overrider = obspar{j};
                
                formula = regexprep(formula, sprintf(rex, overridee), ...
                    num2str(overrider));
            end            
        end
        
        if ~isempty_ext(cvars.keys)
            for j = 1:numel(cvars)
                formula = regexprep(formula, sprintf(rex, cvars.keys(j)), ...
                    num2str(cvars.values{j}));
            end
        end
        
        formula = formula(2:end - 1);
        
        tmp = to_float_if_float(formula);
        if isstring(tmp)
            out.observableFormula{i} = formula;
        else
            out.observableFormula{i} = tmp;
        end
    end
        
    if all(map(@isnumeric, out.observableFormula))
        out.observableFormula = cell2mat(out.observableFormula);
    end
    
    if ismember('observableParameters', out.Properties.VariableNames)
        out.observableParameters = [];
    end
end