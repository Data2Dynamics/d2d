function out = parseOptTable(problem, optTable, p, xdict)
    %TODO
    
    for i = 1:height(optTable)
        x = xdict(string(optTable.time(i)));
        
        obsFormula = parseExpression(problem, ...
            string(optTable.observableFormula{i}));        
        optTable.observableFormula{i} = eval(obsFormula);
        
        if isnumeric(optTable.noiseFormula)
            continue
        end
        
        noiseFormula = parseExpression(problem, ...
            string(optTable.noiseFormula{i}));
        optTable.noiseFormula{i} = eval(noiseFormula);
    end
    
    optTable.observableFormula = cell2mat(optTable.observableFormula);
    mask = strcmp(optTable.Properties.VariableNames, 'observableFormula');
    optTable.Properties.VariableNames{mask} = 'simulation';
    
    if ~isnumeric(optTable.noiseFormula)
        optTable.noiseFormula = cell2mat(optTable.noiseFormula);
    end
    mask = strcmp(optTable.Properties.VariableNames, 'noiseFormula');
    optTable.Properties.VariableNames{mask} = 'sigmas';
    
    out = optTable;
end