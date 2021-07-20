function assert_no_leading_trailing_whitespace(names_list, name)
    %Check that there is no trailing whitespace in elements of an Array or
    %Cell Array of Strings.
    %
    %Arguments:
    %   names_list [string]: 
    %       List of names to check for whitespace.
    %   name string: 
    %       Name of 'names_list' for error messages.
    %
    %Raises:
    %   [TrailingWhitespaceError]
    %       If there is trailing whitespace.
    
    rex = '(?:^\s)|(?:\s$)';
    testfunc = @(s) ~isstring(s) && ~isempty_ext(regexp(s, rex, 'once'));    
    test = map(testfunc, names_list);
    
    if any(test)
        error('%s:TrailingWhitespaceError', name)
    end
end