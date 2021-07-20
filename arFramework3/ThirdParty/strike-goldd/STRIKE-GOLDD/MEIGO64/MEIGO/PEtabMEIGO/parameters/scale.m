function out = scale(parameter, scale_str) %-> numeric
    %Scale parameter according to scale_str.
    %
    %Arguments:
    %   parameter numeric:
    %       Parameter to be scaled.
    %   scale_str string:
    %       One of 'lin' (synonymous with ''), 'log', 'log10'.
    %
    %Returns:
    %   numeric
    %       Parameter scaled according to scale_str.
    
    if ismember(scale_str, ["", "lin"])
        out = parameter;
    elseif strcmp(scale_str, 'log')
        out = log(parameter);
    elseif strcmp(scale_str, 'log10')
        out = log10(parameter);
    else
        error('SCALE:InvalidParScalingError', ...
            'Invalid parameter scaling: %s', scale_str)
    end
end