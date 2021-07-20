function [par_mapping, scale_mapping] = ...
    get_parameter_mapping_for_condition(condition_id, is_preeq, ...
    cur_measurement_df, sbml_model, condition_df, varargin)
    % Create dictionary of parameter value and parameter scale mappings 
    % from PEtab-problem to SBML parameters for the given condition.
    % 
    % Parameters:
    %   condition_id (string, char): 
    %       Condition ID for which to perform mapping.
    %   is_preeq (bool): 
    %       If "True", output parameters will not be mapped.
    %   cur_measurement_df (table): 
    %       Measurement sub-table for current condition
    %   condition_df (table):
    %       PEtab condition DataFrame.
    %   parameter_df *(table/[]):
    %       PEtab parameter DataFrame.
    %   sbml_model (Sbml):
    %       The sbml model with observables and noise specified according to
    %       the PEtab format used to retrieve simulation parameter IDs.     
    %   simulation_parameters *(Dict/[]):
    %       Model simulation parameter IDs mapped to parameter values (output
    %       of "sbml/get_model_parameters(.., with_values=True)").
    %       Optional, saves time if precomputed.
    %   warn_unmapped *(bool/[]):
    %       If "True", log warning regarding unmapped parameters.
    %     
    %   Returns:
    %       ([Dict, Dict])
    %           Tuple of two dictionaries. First dictionary mapping model 
    %           parameter IDs to mapped parameters IDs to be estimated or 
    %           to filled-in values in case of non-estimated parameters.
    %           Second dictionary mapping model parameter IDs to their 
    %           scale.
    %           NaN is used where no mapping exists.
    
    p = inputParser;
    addRequired(p, 'condition_id', @(x) isstring(x) || ischar(x));
    addRequired(p, 'is_preeq', @islogical);
    addRequired(p, 'cur_measurement_df', @istable);
    addRequired(p, 'sbml_model', @Sbml.isSbml);
    addRequired(p, 'condition_df', @istable);
    addParameter(p, 'parameter_df', table.empty(), @istable);
    addParameter(p, 'simulation_parameters', Dict.empty(), @Dict.isDict);
    addParameter(p, 'warn_unmapped', true, @islogical);
    addParameter(p, 'scaled_parameters', false, @islogical);
    
    parse(p, condition_id, is_preeq, cur_measurement_df, sbml_model, ...
        condition_df, varargin{:});
    
    condition_id = p.Results.condition_id;
    is_preeq = p.Results.is_preeq;
    cur_measurement_df = p.Results.cur_measurement_df;
    sbml_model = p.Results.sbml_model;
    condition_df = p.Results.condition_df;
    parameter_df = p.Results.parameter_df;
    simulation_parameters = p.Results.simulation_parameters;
    warn_unmapped = p.Results.warn_unmapped;
    scaled_parameters = p.Results.scaled_parameters;
    
    perform_mapping_checks(cur_measurement_df);
    
    if isempty_ext(simulation_parameters)
        simulation_parameters = get_model_parameters(sbml_model, true);
    end
    
    par_mapping = simulation_parameters.copy;
    
    scale_mapping = Dict(par_mapping.keys, cellstr(repmat("lin", 1, ...
        numel(par_mapping.keys))));  
    
    output_parameters_to_nan(par_mapping);
    
    apply_output_parameter_overrides(par_mapping, cur_measurement_df);
    
    if ~is_preeq
        handle_missing_overrides(par_mapping, warn_unmapped);
    end
    
    apply_condition_parameters(par_mapping, scale_mapping, ...
        condition_id, condition_df, sbml_model);
    
    apply_parameter_table(par_mapping, scale_mapping, parameter_df, ...
        scaled_parameters);    
end