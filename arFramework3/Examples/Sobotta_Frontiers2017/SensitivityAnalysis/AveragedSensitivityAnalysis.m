%% Model specific settings
ar.config.add_c = 200;  % Backward compatibility with old D2D
showLPSA = 1;           % Show the result
calcLPSA = 1;           % Calculate the result
separate = 1;           % Separate figures?
printFig = 0;           % Export the result to a file
socs = 0;               % Show results for SOCS3 output instead of APP
showProfilePoints = 1;  % Show the individual points belonging to the profiles

SensitivityAnalysisOptions;
mkdir('figures/LPSA');

if ( printFig )
    showLPSA = 1;
end

%% Calculate the sensitivities over the profiles
if ( calcLPSA )
    [~,A,B] = intersect(ar.pLabel, pleGlobals.p_labels);
    ar.p(A) = pleGlobals.p(B);
    arSimu(false); arCalcMerit(false); arGetMerit;
    Sref = LPSA( obj, parIDs, 'Perturbation', 0.5, 'ids', 1:length(ids) );
    S=cell(0);
    arPush; arPush;
    bnd = min(cell2mat(cellfun(@min, pleGlobals.chi2s, 'UniformOutput', false))) + pleGlobals.dchi2_point;

    c = 0;
    for jP = 1 : length( pleGlobals.chi2s )
        fprintf( '%d/%d', jP, length(pleGlobals.chi2s) );
        if ~isempty( pleGlobals.ps{jP} )
            c = c + 1;
            first = find(pleGlobals.chi2s{jP} < bnd, 1);
            last = numel(pleGlobals.chi2s{jP}) - find(fliplr(pleGlobals.chi2s{jP}) < bnd, 1) + 1;
            
            ar.p(A) = pleGlobals.ps{jP}(first,B);
            arSimu(false);
            arCalcMerit(false); arGetMerit
            S{c} = LPSA( obj, parIDs, 'Perturbation', 0.5, 'ids', 1:length(ids) );
            c = c + 1;
            ar.p(A) = pleGlobals.ps{jP}(last,B);
            arSimu(false);
            arCalcMerit(false); arGetMerit
            S{c} = LPSA( obj, parIDs, 'Perturbation', 0.5, 'ids', 1:length(ids) );
        end
    end
    arPop;
end

%% Plot results
for a = 1 : length( names )
    appNames{a} = lower(names{a}(1:end-3));
end
if ( isempty( parNames ) )
    parNames = ar.pLabel(parIDs);
end

% LPSA
if ( showLPSA )
    if ( separate )
        for a = 1 : length(names)
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
            [~, p] = setdiff( parIDs, excludePars, 'stable');
            textBar( parNames(p), Sref(p, a) ); hold on;

            xmi = inf;
            xma = 0;
            for C = 1 : length(S)
                plot( mean(S{C}(p,a), 2).', 1:numel(p), 'o' ); hold on;
                xmi = min( [ xmi min( mean(S{C}(p,a),2) ) ] );
                xma = max( [ xma max( mean(S{C}(p,a),2) ) ] );
            end

            title( sprintf( 'Local Parametric Sensitivity Analysis (50%% perturbation) for %s mRNA', appNames{a} ) );
            grid on;

            xlim( [xmi xma] );
            if ( printFig )
                figExport(sprintf('figures/LPSA/LPSA_%s', names{a}));
            end
        end
    else
        % Averaged Sensitivity
        figure;
        cur = lower(names{a}(1:end-3));
        p = 1:numel(parIDs);
        
        % Exclude socs3 since this is not one we wish to suppress
        if ( socs )
            [~,Igenes]=intersect(appNames, 'socs3', 'stable');
        else
            [~,Igenes]=setdiff(appNames, 'socs3', 'stable');
        end
        
        % Determine maxima and minima
        mins = inf * ones(length(p), 1);
        maxs = -inf * ones(length(p), 1);
        for C = 1 : length( S ) 
            mins = min( [ mins, mean(S{C}(p,Igenes),2) ], [], 2 );
            maxs = max( [ maxs, mean(S{C}(p,Igenes),2) ], [], 2 );
        end
        
        if ( showProfilePoints )
            textBar( parNames(p), mean(Sref(p, Igenes), 2), gca, [.2 0 .8] ); hold on;
        else
            textBar( parNames(p), mean(Sref(p, Igenes), 2), gca, [.2 0 .8], mins, maxs ); hold on;
        end
        
        xmi = inf;
        xma = 0;
        for C = 1 : length(S)
            if ( showProfilePoints )
                plot( mean(S{C}(p,Igenes), 2).', 1:numel(p), 'o' ); hold on;
            end
            xmi = min( [ xmi min( mean(S{C}(p,Igenes),2) ) ] );
            xma = max( [ xma max( mean(S{C}(p,Igenes),2) ) ] );
        end
        
        title( sprintf( 'Averaged Local Parametric Sensitivity Analysis (50%% perturbation) for APP mRNA expression' ) );
        grid on;
        
        dff = abs(xma - xmi);
        xlim( [xmi-0.1*dff xma+0.1*dff] );
        
        if ( printFig )
            if ( socs )
                figExport('figures/LPSA/LPSA_Averaged_SOCS3');
            else
                figExport('figures/LPSA/LPSA_Averaged_APP');
            end
        end

    end
end