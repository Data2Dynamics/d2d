function llh = calculateSingleLlh(measurement, simulation, scale, ...
    noise_distribution, noise_value)
    %Calculates single log-likelihood.
    %
    %Arguments:
    %   measurement numeric:
    %       Observable measurement unscaled.
    %   simulation numeric:
    %       Observable simulation unscaled.
    %   scale string:
    %       Observable scale, "lin", "log" or "log10".
    %   noise_distribution string:
    %       Observable noise distribution, "normal" or "laplace".
    %   noise_value numeric:
    %       Observable noise value.
    %
    %Returns:
    %   numeric
    %       Single llh.
    
    m = measurement;
    s = simulation;
    sigma = noise_value;
    
    if strcmp(noise_distribution, 'normal') && strcmp(scale, 'lin')
        nllh = 0.5*log(2*pi*sigma^2) + 0.5*((s-m)/sigma)^2;
    elseif strcmp(noise_distribution, 'normal') && strcmp(scale, 'log')
        nllh = 0.5*log(2*pi*sigma^2*m^2) + 0.5*((log(s)-log(m))/sigma)^2;
    elseif strcmp(noise_distribution, 'normal') && strcmp(scale, 'log10')
        nllh = 0.5*log(2*pi*sigma^2*m^2*log(10)^2) + ...
            0.5*((log10(s)-log10(m))/sigma)^2;
    elseif strcmp(noise_distribution, 'laplace') && strcmp(scale, 'lin')
        nllh = log(2*sigma) + abs((s-m)/sigma);
    elseif strcmp(noise_distribution, 'laplace') && strcmp(scale, 'log')
        nllh = log(2*sigma*m) + abs((log(s)-log(m))/sigma);
    elseif strcmp(noise_distribution, 'laplace') && strcmp(scale, 'log10')
        nllh = log(2*sigma*m*log(10)) + abs((log10(s)-log10(m))/sigma);
    end
    
    llh = - nllh;
end