function out = calculateSingleChi2(measurement, simulation, noise_value, scale)
    %Calculates single chi2 as:
    %
    %          chi2_i = (simulation_i - measure_i)^2 / noise_value
    %
    %Arguments:
    %   measurement numeric:
    %       Observable measurement unscaled.
    %   simulation numeric:
    %       Observable simulation unscaled.
    %   noise_value numeric:
    %       Observable noise value.
    %   scale string:
    %       Observable scale, "lin", "log" or "log10".
    %
    %Returns:
    %   numeric
    %       Single chi2.
    
    if strcmp(scale, 'lin')
        out = ((simulation - measurement)/noise_value)^2;
    elseif strcmp(scale, 'log')
        out = ((log(simulation) - log(measurement))/noise_value)^2;
    elseif strcmp(scale, 'log10')
        out = ((log10(simulation) - log10(measurement))/noise_value)^2;
    end
end