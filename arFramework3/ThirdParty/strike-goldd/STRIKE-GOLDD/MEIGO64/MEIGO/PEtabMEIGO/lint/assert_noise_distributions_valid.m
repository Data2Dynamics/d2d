function assert_noise_distributions_valid(observables_df)
    %Ensure that noise distributions and transformation for observables are
    %valid.
    %
    %Arguments:
    %   observables_df [table]:
    %       PEtab observables table.
    %
    %Raises:
    %    [ObsTrafoInvalidError, NoiseDistrInvalidError]
    %       If noise distribution or observable transformation is invalid.
    
    OBS_TRAFOS = ["", "lin", "log", "log10"];
    NOISE_DISTR = ["", "normal", "laplace"];  
    
    cols = string(observable.Properties.VariableNames);    
    if ismember('observableTransformation', cols)
        trafo = observables_df.observableTransformation;
        testfunc = @(x) ~(isnumeric(x) || ismember(x, OBS_TRAFOS));
        
        test = map(testfunc, trafo);
        if any(test)
            error('TRAFO:ObsTrafoInvalidError', ['Unrecognized ' ...
                'observable transformation in observable table'])
        end
    end
    
    if ismember('noiseDistribution', cols)
        noise = observables_df.noiseDistribution;
        testfunc = @(x) ~(isnumeric(x) || ismember(x, NOISE_DISTR));
        
        test = map(testfunc, noise);
        if any(test)
            error('DISTR:NoiseDistrInvalidError', ['Unrecognized ' ...
                'noise distribution in observable table'])
        end
    end
end