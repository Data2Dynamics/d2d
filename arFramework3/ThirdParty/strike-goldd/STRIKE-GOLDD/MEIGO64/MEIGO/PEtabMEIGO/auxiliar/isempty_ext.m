function out = isempty_ext(input) %-> bool
    %Returns true if "input" is empty. Matlab builtin isempty extension
    %which considers character array '' as an empty object. 
    
    out = isempty(input) || isequal(input, '');
end