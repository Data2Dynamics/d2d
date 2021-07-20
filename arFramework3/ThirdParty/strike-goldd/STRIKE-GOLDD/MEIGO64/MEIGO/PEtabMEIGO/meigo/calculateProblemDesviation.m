function [llh, chi2, simtable] = calculateProblemDesviation(problem, p, ...
    sims)
    %Calculates desviation measures (chi2 and log-likelihood) for a problem
    %described in PEtab format.
    %
    %Arguments:
    %   problem Problem:
    %       Problem class object.
    %   p [numeric]:
    %       Vector of free parameters.
    %   sims Optional bool/*false:
    %       Whether to return simulations table.
    %
    %
    
    
    if nargin == 2
        sims = false;
    end
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    
    if ~sims
        simtable = [];
    end
    
    chi2 = 0;
    llh = 0;
    
    status = writeDynFunc(problem, problem.folder, problem.sbml_model.id);
    
    simcond = problem.get_simulation_conditions_from_measurement_df;
    simcond.Count = [];
    
    is_preeq = ismember('preequilibrationConditionId', ...
        simcond.Properties.VariableNames);
    
    if sims
        if is_preeq
            simulations = problem.measurement_df(:, ["observableId" ...
                "preequilibrationConditionId" "simulationConditionId" ...
                "time"]);
        else
            simulations = problem.measurement_df(:, ["observableId" ...
                "simulationConditionId" "time"]);
        end
        
        simtable = simulations;
        simtable.simulation = zeros(height(simtable), 1);
        simstruct = table2struct(simulations);
    end
    
    x0 = Dict();
    simids = string(simcond.simulationConditionId);
    if is_preeq
        fprintf('Calculating preequilibration...\n')
        
        preeqids = string(unique(simcond.preequilibrationConditionId));        
        
        for i = 1:numel(preeqids)
            c = cvarsForCondition(problem, preeqids(i));
            x0i = eval(str2sym(parsex0(problem, preeqids(i))));
            
            [~, x0(preeqids(i))] = simulateForCondition(problem, ...
                preeqids(i), x0i, p, true);
        end
        
        x0i = string(simcond.preequilibrationConditionId);
        
        fprintf('Preequilibration succesfully completed.\n')
    else     
        for i = 1:numel(simids)
            c = cvarsForCondition(problem, simids(i));
            x0(simids(i)) = eval(str2sym(parsex0(problem, simids(i))));
        end
        
        x0i = string(simcond.simulationConditionId);
    end
    
    fprintf('Simulating...\n')
    for i = 1:height(simcond)
        if is_preeq
            for j = 1:numel(x0i)
                x0(x0i(j)) = eval(str2sym(parsex0(problem, simids(i), ...
                    string(x0(x0i(j))))));
            end
        end
        
        [t, x] = simulateForCondition(problem, simids(i), x0(x0i(i)), p);     
        
        xdict = Dict();
        for j = 1:numel(t)
            xdict(string(t(j))) = x(j, :);
        end
        
        optTable = getOptimizationTableForCondition(problem, simids(i));
        optTable = parseOptTable(problem, optTable, p, xdict);
        
        for j = 1:height(optTable)
            llh = llh + calculateSingleLlh(optTable.measurement(j), ...
                optTable.simulation(j), ...
                optTable.observableTransformation{j}, ...
                optTable.noiseDistribution{j}, optTable.sigmas(j));
            
            chi2 = chi2 + calculateSingleChi2(optTable.measurement(j), ...
                optTable.simulation(j), optTable.sigmas(j), ...
                optTable.observableTransformation{j});
        end        
        
        if sims
            uniqcond = problem.get_simulation_conditions_from_measurement_df;
            for j = 1:height(optTable)
                row = struct();
                
                row.observableId = optTable.observableId{j};
                row.simulationConditionId = ...
                    uniqcond.simulationConditionId{i};
                row.time = optTable.time(j);
                
                if is_preeq
                    row.preequilibrationConditionId = ...
                        uniqcond.preequilibrationConditionId{i};
                end
                
                mask = false(1, numel(simstruct));
                for k = 1:numel(simstruct)
                    mask(k) = isequal(simstruct(k), row);
                end
                
                simtable.simulation(mask) = optTable.simulation(j);
            end
        end
    end
    
    fprintf('Simulation succesfully completed.\n')
    
    delete(fullfile(problem.folder, [problem.sbml_model.id '_dyn.m']));
end