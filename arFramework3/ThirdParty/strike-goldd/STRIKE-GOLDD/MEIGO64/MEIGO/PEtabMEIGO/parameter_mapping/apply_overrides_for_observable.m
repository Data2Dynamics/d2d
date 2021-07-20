function apply_overrides_for_observable(mapping, observable_id, ...
    override_type, overrides)
    %Apply parameter-overrides for observables and noises to mapping
    %matrix.
    %
    %Parameters:
    %	mapping (Dict): 
    %   	mapping dict to which to apply overrides.
    %	observable_id (string, char): 
    %   	observable ID.
    %   override_type (string, char): 
    %   	'observable' or 'noise'
    %	overrides ([string]): 
    %   	list of overrides for noise or observable parameters.
    
    for i = 1:numel(overrides)
        override = overrides{i};
        
        overridee_id = sprintf('%sParameter%d_%s', override_type, i, ...
            observable_id);
        
        try
            mapping(overridee_id) = override;
        catch
            error('APPLY_OVERRIDES_FOR_OBSERVABLE:TypeError', ...
                ['Cannot override %s parameter %s for observable %s.' ...
                'Ensure there exists an %d definition containing the ' ...
                'correct number of placeholder parameters.'], ...
                override_type, overridee_id, observable_id, override_type)
        end
    end
end