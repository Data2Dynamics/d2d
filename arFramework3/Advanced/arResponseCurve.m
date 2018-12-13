% rate = arResponseCurve( name, indep1, indep2, varargin )
%
% Show instantaneous response curve w.r.t. two variables of a model.
%
%   name 			Name of the flux.
%	indep1          Independent variable 1 (x-axis).
% 	indep2          Independent variable 2.
%   range1          Specify range for variable 1 manually.
%   range2          Specify range for variable 2 manually.
%
% This function can be used to probe how a specific flux expression depends
% on two variables. The flux is evaluated for specific values for the
% variables specified by indep1 and indep2 over the ranges range1 and
% range2. For all the variables that are not specified, state and parameter
% values are taken from a reference condition (1 by default, but can be
% changed by optional arguments) at a reference time point (initial by 
% default). This function is primarily intended for enzymes, where it can
% show how the rate depends on the specific states or parameters.
%
%   Optional arguments:
%   'timepoints'    Specify simulation time point to use for values of the 
%                   other states (default is initial condition).
%	'model'         Followed by value which specifies which model to use 
%                   for reference species.
%	'condition'     Followed by index which indicates the condition that
%                   should be used to obtain reference values for 
%                   species that are not sampled.
%	'relative'      Divide rate by rate which is obtained by taking
%                   the limit of the substrates to infinity and 
%                   products to 0 (for enzymes this'd typically be vmax)
%
% See also arFractionalFlux

function rate = arResponseCurve( name, indep1, indep2, varargin )

    global ar;
    warning('off', 'symbolic:sym:sym:DeprecateExpressions');
    ylog = 0;
    m = 1;
    cond = 1;
    timepoints = 1;
    miniTresh = 1e-16;
    mRange = 5;
    
    args = {'timepoints', 'condition', 'model', 'relative', 'range', 'noclear', 'hold', 'custom', 'range1', 'range2'};
    extraArgs = [1, 1, 1, 0, 1, 0, 0, 1, 1, 1];
    opts = argSwitch( args, extraArgs, {}, 0, varargin );

    if opts.timepoints
        timepoints = opts.timepoints_args;
    end
    if opts.model
        model = opts.model_args;
        if ( ~isnumeric( model ) || ( model > numel( ar.model ) ) )
            error( 'Invalid argument passed for model' );
        end
    end
    if opts.condition
        cond = opts.condition_args;
        if ( ~isnumeric( cond ) || ( cond > numel( ar.model(model).condition ) ) )
            error( 'Invalid argument passed for condition' );
        end    
    end
    if opts.range
        mRange = opts.range_args;
    end
    
    arSubplot(1,1,1); 
    if ( (~opts.noclear) && (~opts.hold) )
        cla;
    end
    hold on;
    NX = ceil(sqrt(numel(timepoints)));
    NY = ceil(numel(timepoints)/NX);
    try
        arSimu(false, true, true);
    catch
        disp( 'Simulation failure' );
    end
    for ti = 1 : numel( timepoints )
        
        if ( numel( timepoints ) > 1 )
            subplot(NX, NY, ti); hold on;
        end
        tp = timepoints(ti);
        
        enzyme = ismember( ar.model(m).v, name );
        
        fvsym = sym( ar.model(m).fv{enzyme} );
        func = subs( fvsym, ar.model(m).condition(cond).pold, ar.model(m).condition(cond).fp.' );
        func = simplify( func );
        
        pLabels = ar.pLabel;
        xLabels = ar.model(m).x;
        zLabels = ar.model(m).z;
        uLabels = ar.model(m).u;
        
        pValues = arGetPars( ar.pLabel, 0 );
        xValues = ar.model(m).condition(cond).xFineSimu(tp, :) + 0;
        zValues = ar.model(m).condition(cond).zFineSimu(tp, :) + 0;
        uValues = ar.model(m).condition(cond).uFineSimu(tp, :) + 0;
        
        labels = [ pLabels, xLabels, zLabels, uLabels ];
        values = [ pValues, xValues, zValues, uValues ];
        
        [~, ~, Iref] = intersect( {indep1, indep2}, labels, 'stable' );
        refValues = values(Iref);
        
        % Overrides
        if ( opts.custom )
            for a = 1 : 2 : numel( opts.custom_args )
                values( strcmp( labels, opts.custom_args{a} ) ) = opts.custom_args{a+1};
            end
        end
        
        % Remove the variables we want to scan
        [labels, I] = setdiff( labels, {indep1, indep2} );
        values = values(I);

        if ( opts.relative )
            % Get the vmax of this reaction
            substrates = ar.model(m).fv_source{enzyme};
            products = ar.model(m).fv_target{enzyme};
            products = setdiff( products, substrates );
            limF = func;
            for i = 1 : numel( substrates )
                limF = limit( limF, substrates{i}, inf );
            end
            for i = 1 : numel( products )
                limF = limit( limF, products{i}, 0 );
            end
            
            disp( 'Limit value (verify whether this contains only vmax expression):' );
            
            % Normalize by vmax
            func = func / limF;
            
            rateName = sprintf( 'v_{%s} / v_{%s}^{max}', ar.model(m).v{enzyme}, ar.model(m).v{enzyme} );
        else
            rateName = sprintf( 'v_{%s}', ar.model(m).v{enzyme} );
        end
        
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
        if ( opts.range1 )
            responseRange1 = opts.range1_args;
        else
            responseRange1 = 10.^[-mRange : .05 : mRange]; %#ok
        end
        if ( opts.range2 )
            responseRange2 = opts.range2_args;
        else
            responseRange2 = 10.^[-mRange : .5  : mRange]; %#ok
        end

        rate = zeros( numel( responseRange2 ), numel( responseRange1 ) );
        cmap = parula( numel( responseRange2 ) );

        doses = zeros( numel(responseRange2 ) + 1, 1 );
        plots = zeros( numel(responseRange2 ) + 1, 1 );
        for c1 = 1 : numel( responseRange2 )
            rate(c1,:) = mFunc( responseRange1, responseRange2(c1) );
            if ( ylog )
                rate(c1,rate(c1,:)<miniTresh) = miniTresh;
            end
            doses(c1) = responseRange2(c1);
            plots(c1) = plot( responseRange1, rate(c1, :), 'Color', cmap(c1, :) );
        end
        r = mFunc( responseRange1, refValues(2) );
        r(:, r<miniTresh) = miniTresh;
        doses(end) = refValues(2);
        plots(end) = plot( responseRange1, r, 'k', 'LineWidth', 2 );
        plot( [refValues(1), refValues(1)], [min(min(rate)), max(max(rate))], 'k--' );

        xlabel( strrep( indep1, '_', '\_' ) );
        ylabel( sprintf( '%s flux [%s]', enzyme, ar.model(m).vUnits{enzyme,2} ) );
        if ( ylog )
            ylim( [max([1e-12, min(min(rate))]), max(max(rate))] );
        else
            ylim( [min(min(rate)), max(max(rate))] );
        end

        set( gca, 'XScale', 'log' )
        if ( ylog )
            set( gca, 'YScale', 'log' )
            ylabel( sprintf( 'log_{10}(%s)', rateName ) );
        else
            ylabel( sprintf( '%s', rateName ) );
        end
        
        set( gca, 'CLim', [min(log10(responseRange2)), max(log10(responseRange2))] );
        c = colorbar;
        colormap(parula);
        ylabel( c, sprintf( 'log_{10}(%s)', strrep(indep2, '_', '\_') ) );
    
        % If the interactivity system is enabled, register the callbacks
        % and provide arInteractivity with the required data.
        if ( exist( 'arInteractivity', 'file' ) )
            if ( arInteractivity )
                arInteractivity( 'arResponseCurve', plots, doses, indep2 );
            end
        end        
        
        if ( nargout == 0 )
            clear rate;
        end
    end
    
    warning('on', 'symbolic:sym:sym:DeprecateExpressions');
end
