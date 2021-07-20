function [preeq, sim] = get_optimization_to_simulation_parameter_mapping( ...
    condition_df, measurement_df, varargin) %-> [[Dict], [Dict]]
    %Create list of mapping dicts from PEtab-problem to SBML
    %parameters.
    %
    %Mapping can be performed in parallel. The number of threads is controlled
    %by the environment variable with the name of petab.ENV_NUM_THREADS.
    %
    %Parameters:
    %   condition_df, measurement_df table:
    %       The dataframes in the PEtab format.
    %   parameter_df Optional table/*[]:
    %       Parameter table in PEtab format.
    %   observable_df Optional table/*[]:
    %       Observables table in PEtab format.
    %   sbml_model Optional libsbml struct/*[]:
    %       The sbml model with observables and noise specified
    %       according to the PEtab format.
    %   simulation_conditions Optional table/*[]:
    %       Table of simulation conditions as created by
    %       "petab.get_simulation_conditions".
    %   warn_unmapped Optional bool/*true:
    %       If "True", log warning regarding unmapped parameters
    %   scaled_parameters Optional bool/*false:
    %       Whether parameter values should be scaled.
    %
    %Returns:
    %   [[Dict], [Dict]]:
    %       Parameter value and parameter scale mapping for all
    %       conditions.
    
    preeq = Dict();
    sim = Dict();
    
    p = inputParser;
    addRequired(p, 'condition_df', @istable);
    addRequired(p, 'measurement_df', @istable);
    addParameter(p, 'parameter_df', table.empty(), @istable);
    addParameter(p, 'observable_df', table.empty(), @istable)
    addParameter(p, 'sbml_model', Sbml.empty(), @Sbml.isSbml)
    addParameter(p, 'simulation_conditions', table.empty(), @istable)
    addParameter(p, 'warn_unmapped', true, @islogical)
    addParameter(p, 'scaled_parameters', false, @islogical)
    parse(p, condition_df, measurement_df, varargin{:});
    
    condition_df = p.Results.condition_df;
    measurement_df = p.Results.measurement_df;
    parameter_df = p.Results.parameter_df;
    observable_df = p.Results.observable_df;
    sbml_model = p.Results.sbml_model;
    simulation_conditions = p.Results.simulation_conditions;
    warn_unmapped = p.Results.warn_unmapped;
    scaled_parameters = p.Results.scaled_parameters;
    
    perform_mapping_checks(measurement_df)
    
    if isempty_ext(simulation_conditions)
        simulation_conditions = get_simulation_conditions(measurement_df);
    end
    
    simulation_parameters = get_model_parameters(sbml_model, true);    
    
    if ~isempty_ext(observable_df)
        output_parameters = get_output_parameters(observable_df, ...
            sbml_model);
        
        simulation_parameters.addpairs(string(output_parameters), ...
            NaN(1, numel(output_parameters)));
    end
    
    tmp = cell(height(simulation_conditions), 4);
    for i = 1:height(simulation_conditions)
        condition = simulation_conditions(i, :);
        
        [tmp{i, :}] = map_condition(condition, measurement_df, ...
            condition_df, parameter_df, sbml_model, ...
            simulation_parameters, warn_unmapped, scaled_parameters);
    end
    
    if ismember('preequilibrationConditionId', ...
            simulation_conditions.Properties.VariableNames)
        
        for i = 1:height(simulation_conditions)
            condid = simulation_conditions.preequilibrationConditionId{i};            
            preeq(condid) = tmp(i, 1:2);
        end
    end
    
    for i = 1:height(simulation_conditions)
        condid = simulation_conditions.simulationConditionId{i};
        sim(condid) = tmp(i, 3:4);
    end
end
