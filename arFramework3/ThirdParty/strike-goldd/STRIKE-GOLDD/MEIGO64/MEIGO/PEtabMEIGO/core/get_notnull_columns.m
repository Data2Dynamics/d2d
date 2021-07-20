function out = get_notnull_columns(df, candidates) %-> [string]
    %Return list of 'df' columns in 'candidates' which are not all NaN.
    %
    %Parameters:
    %   df [table]:
    %       MATLAB table.
    %   candidates [string]:
    %       Columns of 'df' to consider.
    %
    %Returns:
    %   [string]
    %       String array of columns of 'df' in 'candidates' which are not
    %       all NaN.
    
    df_cols = string(df.Properties.VariableNames);
    [candidates, ia] = intersect(df_cols, candidates, 'stable');
    
    testfunc = @(c) ~(isnumeric(c) && all(isnan(c)));    
    out = candidates(varfun(testfunc, df(:, ia), 'OutputFormat', ...
        'uniform'));
end