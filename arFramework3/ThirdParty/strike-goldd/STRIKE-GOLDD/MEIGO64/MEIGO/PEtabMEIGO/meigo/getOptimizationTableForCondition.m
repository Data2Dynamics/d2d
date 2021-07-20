function out = getOptimizationTableForCondition(problem, condid)
    %TODO
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    observable_df = problem.observable_df;
    observable_df.Properties.RowNames = observable_df.observableId;
    observable_df.observableId = [];
    
    mask = strcmp(problem.measurement_df.simulationConditionId, condid);    
    out = problem.measurement_df(mask, ["observableId" ...
        "measurement" "time"]);
    
    out.observableFormula = map( ...
        @(o)observable_df.observableFormula(o), out.observableId);
    
    out.noiseFormula = map(@(o)observable_df.noiseFormula(o), ...
        out.observableId);
         
    if ismember('observableParameters', ...
            problem.measurement_df.Properties.VariableNames)
        
        out.observableParameters = ...
            problem.measurement_df.observableParameters(mask);
    end
    
    if ismember('noiseParameters', ...
            problem.measurement_df.Properties.VariableNames)
        
        out.noiseParameters = ...
            problem.measurement_df.noiseParameters(mask);
    end
         
    if ismember('observableTransformation', ...
            problem.observable_df.Properties.VariableNames)
        
        out.observableTransformation = map( ...
            @(o)observable_df.observableTransformation(o), ...
            out.observableId);
    else
        out.observableTransformation = ...
            cellstr(repmat("lin", height(out), 1));
    end
    
    if ismember('noiseDistribution', ...
            problem.observable_df.Properties.VariableNames)
        
        out.noiseDistribution = map( ...
            @(o)observable_df.noiseDistribution(o), out.observableId);
    else
        out.noiseDistribution = ...
            cellstr(repmat("normal", height(out), 1));        
    end
    
    out = parseObservableFormulas(problem, out, condid);
    out = parseNoiseFormulas(problem, out, condid);    
end