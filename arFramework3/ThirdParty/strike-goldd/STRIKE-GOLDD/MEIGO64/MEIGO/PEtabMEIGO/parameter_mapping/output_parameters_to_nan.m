function output_parameters_to_nan(mapping)
    %Set output parameters in mapping dictionary to NaN.
    
    rex = '^(noise|observable)Parameter[0-9]+_';
    mask = map(@(s) ~isempty_ext(regexp(s, rex)), mapping.keys);
    
    tmp = mapping.keys(mask);
    mapping.addpairs(tmp, NaN(1, numel(tmp)));
end