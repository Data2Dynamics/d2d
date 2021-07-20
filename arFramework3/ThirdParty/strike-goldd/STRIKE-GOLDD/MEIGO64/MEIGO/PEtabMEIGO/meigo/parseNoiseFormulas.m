function out = parseNoiseFormulas(problem, opt_problem, condid)
    %TODO   
    
    out = opt_problem;
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    if ~ismember('noiseParameters', out.Properties.VariableNames)
        return;
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
        formula = [' ' out.noiseFormula{i} ' '];
        
        if ismember('noiseParameters', out.Properties.VariableNames)
            if isnumeric(out.noiseParameters)
                noisepar = out.noiseParameters(i);
            else
                noisepar = out.noiseParameters{i};
            end
            noisepar = split_parameter_replacement_list(noisepar);            

            for j = 1:numel(noisepar)
                overridee = sprintf('noiseParameter%d_%s', j, obsid);
                overrider = noisepar{j};
                
                formula = regexprep(formula, sprintf(rex, overridee), ...
                    num2str(overrider));
            end            
        end
        
        if ~isempty(cvars.keys)
            for j = 1:numel(cvars)
                formula = regexprep(formula, sprintf(rex, cvars.keys(j)), ...
                    num2str(cvars.values{j}));
            end
        end
        
        formula = formula(2:end - 1);
        
        tmp = to_float_if_float(formula);
        if isstring(tmp)
            out.noiseFormula{i} = formula;
        else
            out.noiseFormula{i} = tmp;
        end
    end
    
    if all(map(@isnumeric, out.noiseFormula))
        out.noiseFormula = cell2mat(out.noiseFormula);
    end
    
    if ismember('noiseParameters', out.Properties.VariableNames)
        out.noiseParameters = [];
    end
end