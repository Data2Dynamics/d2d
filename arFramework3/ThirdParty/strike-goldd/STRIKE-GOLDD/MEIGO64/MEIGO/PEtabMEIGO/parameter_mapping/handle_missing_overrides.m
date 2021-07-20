function handle_missing_overrides(mapping_par_opt_to_par_sim, warn, ...
    condition_id)
    %Find all observable parameters and noise parameters that were not 
    %mapped and set their mapping to NaN.
    %
    %Assumes that parameters matching "(noise|observable)Parameter[0-9]+_" 
    %were all supposed to be overwritten.
    %
    %Parameters:
    %    mapping_par_opt_to_par_sim (Dict):
    %        Output of get_parameter_mapping_for_condition
    %    warn *(bool/true):
    %        If True, log warning regarding unmapped parameters
    %    condition_id *((string, char)/[]):
    %        Optional condition ID for more informative output
    
    if nargin == 2
        condition_id = string.empty();
    elseif nargin == 1
        warn = true;
        condition_id = string.empty();
    end
    
    rex = '^(noise|observable)Parameter[0-9]+_';
    
    keys = mapping_par_opt_to_par_sim.keys;
    vals = mapping_par_opt_to_par_sim.values;
    
    mask = map(@(s) ~isnumeric(s) && ~isempty_ext(regexp(s, rex)), vals);
    
    missed_vals = keys(mask);
    mapping_par_opt_to_par_sim.addpairs(keys(mask), vals(mask));
    
    if ~isempty(missed_vals) && warn
        warning('HANDLE_MISSING_OVERRIDES:UnmappedParsWarning', ...
            ['Could not map all overrides for condition %s. Usually, ' ...
            'this is just due to missing data points'], condition_id);
    end
end