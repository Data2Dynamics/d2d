%% Model specific settings
ar.config.add_c = 200;
calcLPSA = 1;

SensitivityAnalysisOptions;

%% Calculate the sensitivities
if ( calcLPSA )
    S = LPSA( obj, parIDs, 'Perturbation', 0.5, 'ids', 1:length(ids) );
end

%% Plot results
for a = 1 : length( names )
    appNames{a} = lower(names{a}(1:end-3));
end
if ( isempty( parNames ) )
    parNames = ar.pLabel(parIDs);
end

% LPSA
if ( calcLPSA )
    for a = 1 : length( names )
        figure;
        % Exclude all APP genes except the current one
        cur = lower(names{a}(1:end-3));
        other = setdiff(setdiff(appNames, cur), 'socs3');

        excludePars = [];
        for c = 1 : length( other )
            excludePars = union( excludePars, arFindPar([other{c} 'rna_delay']) );
            excludePars = union( excludePars, arFindPar([other{c} 'rna_hill']) );
            excludePars = union( excludePars, arFindPar([other{c} 'rna_des']) );
            excludePars = union( excludePars, arFindPar([other{c} 'rna_ka']) );
            excludePars = union( excludePars, arFindPar([other{c} 'rna_pro']) );
        end
        [~, p] = setdiff( parIDs, excludePars );
        textBar( parNames(p), S(p, a) );
        title( sprintf( 'Local Parametric Sensitivity Analysis (50%% perturbation) for %s mRNA', appNames{a} ) );
        grid on;
    end
end