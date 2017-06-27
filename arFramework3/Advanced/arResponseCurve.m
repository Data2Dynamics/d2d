% Show rate response curve w.r.t. two variables of a model.
%
% Usage:
%		arResponseCurve( name, indep1, indep2, timepoints )
%
% Arguments:
% 		name 			- Name of the flux
%		indep1		- Independent variable 1 (x-axis)
% 		indep2 		- Independent variable 2
%		timepoints 	- Simulation time point to use for values of the other states (default is initial condition)
%
% NOTE: THIS FUNCTION IS CURRENTLY IN ALPHA STATUS!
% PLEASE DO NOT EDIT THIS FUNCTION YET. IT IS STILL UNDER ACTIVE DEVELOPMENT.
% FUNCTION ARGUMENTS ARE LIKELY SUBJECT TO CHANGE IN THE NEAR FUTURE.
function rate = arResponseCurve( name, indep1, indep2, timepoints )

    global ar;
    
    ylog = 0;
    m = 1;
    cond = 1;
    miniTresh = 1e-16;
    if ~exist( 'timepoints', 'var' )
        timepoints = 1;
    end
    
    figure; hold on;
    NX = ceil(sqrt(numel(timepoints)));
    NY = ceil(numel(timepoints)/NX);
    for ti = 1 : numel( timepoints )
        
        subplot(NX, NY, ti); hold on;
        tp = timepoints(ti);
        
        enzyme = ismember( ar.model.v, name );
        func = sym( ar.model.fv{enzyme} );

        pLabels = ar.pLabel;
        xLabels = ar.model.x;
        zLabels = ar.model.z;
        uLabels = ar.model.u;
        pValues = arGetPars( ar.pLabel, 0 );
        xValues = ar.model(m).condition(cond).xFineSimu(tp, :);
        zValues = ar.model(m).condition(cond).zFineSimu(tp, :);
        uValues = ar.model(m).condition(cond).uFineSimu(tp, :);
        
        labels = [ pLabels, xLabels, zLabels, uLabels ];
        values = [ pValues, xValues, zValues, uValues ];      

        [~, ~, Iref] = intersect( labels, {indep1, indep2} );
        refValues = values(Iref);

        % Remove the variables we want to scan
        [labels, I] = setdiff( labels, {indep1, indep2} );
        values = values(I);

        % Substitute variables
        func = subs( func, labels, values );

        vars = symvar( func );
        if ( sum( ismember( vars, indep1 ) ) == 0 )
            error( 'Independent variable %s not found in rate equation %s', indep1, char(func) );
        end
        if ( sum( ismember( vars, indep2 ) ) == 0 )
            error( 'Independent variable %s not found in rate equation %s', indep2, char(func) );
        end

        % Generate the matlab function
        mFunc = matlabFunction( func, 'vars', {indep1, indep2} );

        % Get reference amounts
        responseRange1 = 10.^[-5 : .1 : 5];
        responseRange2 = 10.^[-5 : .5 : 5];

        rate = zeros( numel( responseRange2 ), numel( responseRange1 ) );
        cmap = parula( numel( responseRange2 ) );

        for c1 = 1 : numel( responseRange2 )
            rate(c1,:) = mFunc( responseRange1, responseRange2(c1) );
            if ( ylog )
                rate(c1,rate(c1,:)<miniTresh) = miniTresh;
            end
            plot( responseRange1, rate(c1, :), 'Color', cmap(c1, :) );
        end
        r = mFunc( responseRange1, refValues(2) );
        r(:, r<miniTresh) = miniTresh;
        plot( responseRange1, r, 'k', 'LineWidth', 2 );
        plot( [refValues(1), refValues(1)], [min(min(rate)), max(max(rate))], 'k--' );

        xlabel( indep1 );
        ylabel( sprintf( '%s flux [%s]', enzyme, ar.model.vUnits{enzyme,2} ) );
        if ( ylog )
            ylim( [max([1e-12, min(min(rate))]), max(max(rate))] );
        else
            ylim( [min(min(rate)), max(max(rate))] );
        end

        set( gca, 'XScale', 'log' )
        if ( ylog )
            set( gca, 'YScale', 'log' )
        end
        set( gca, 'CLim', [min(log10(responseRange2)), max(log10(responseRange2))] );
        c = colorbar;
        colormap(parula);
        ylabel( c, sprintf( 'log_{10}(%s)', indep2 ) );
    
    end
end
