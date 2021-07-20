function perform_mapping_checks(measurement_df)
    %Check for PEtab features which we can't account for during parameter
    %mapping.
    
    if measurement_table_has_timepoint_specific_mappings(measurement_df)
        error('MAPPINGCHECKS:ValueError', ['Timepoint-specific ' ...
            'parameter overrides currently unsupported']);
    end
end