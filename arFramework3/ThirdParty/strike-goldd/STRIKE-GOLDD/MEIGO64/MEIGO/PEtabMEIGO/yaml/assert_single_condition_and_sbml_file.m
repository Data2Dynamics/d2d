function assert_single_condition_and_sbml_file(problem_config)
    %Check that there is only a single condition file and a single SBML
    %file specified.
    %
    %Arguments:
    %   problem_config [struct]:
    %       Structure as defined in the YAML schema inside the `problems`
    %       list.
    % Raises:
    %   [NotImplementedError]
    %       If multiple condition or SBML files specified.
    
    if numel(problem_config.sbml_files) > 1 || ...
            numel(problem_config.condition_files) > 1
        
        error('YAML:NotImplementedError', ['Support for multiple ' ...
            'models or condition files is not yet implemented'])
    end
end