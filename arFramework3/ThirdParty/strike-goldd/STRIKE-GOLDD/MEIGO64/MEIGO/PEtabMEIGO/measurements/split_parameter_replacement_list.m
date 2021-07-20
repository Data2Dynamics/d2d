function out = split_parameter_replacement_list(list_string, delim) %-> string
    %Split values in observableParameters and noiseParameters in
    %measurement table.
    %
    %Arguments:
    %   list_string string:
    %       delim-separated stringified list.
    %   delim Optional string/*[";"]:
    %       delimiter.
    %
    %Returns:
    %   cell
    %       List of split values. Numeric values converted to float.    
    
    if nargin == 1
        delim = ';';
    end
    
    if isnumeric(list_string)
        if isnan(list_string)
            out = [];
        else        
            out = {list_string};
        end
    else
        tmp = split(list_string, delim);
        out = transpose(map(@to_float_if_float, tmp, 'UniformOutput', ...
            false));
    end
end