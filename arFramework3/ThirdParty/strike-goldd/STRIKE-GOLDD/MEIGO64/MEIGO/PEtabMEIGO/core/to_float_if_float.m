function out = to_float_if_float(x) %-> any
    %Return input as float if possible, otherwise return as is.
    %
    %Parameters:
    %   x [any]:
    %       Anything.
    %
    %Returns:
    %   "x" as float if possible, otherwise "x".
    
    if iscell(x)
        x = x{:};
    end
    
    if ischar(x)
        x = string(x);
    end
    
    out = double(x);
    if isstring(x) && isnan(out)
        out = x;
    end            
end