function h1 = qPlotDR( ar, dataSets, variable, tp, color, state, infinityLocation, selectOne )
    
    global errorBarMode;

    opts = {'Color', color, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', .5};    
    optsNoMarker = {'Color', min(color + 0.4, [1,1,1])};

    if ( nargin < 8 )
        selectOne = [];
    end
       
    if ( nargin < 6 ) || ( isempty(state) )
        if (length( ar.model.data( dataSets(1) ).y ) > 1)
            error( 'qPlotDR Error: Specify state' );
        else
            state = ar.model.data( dataSets(1) ).y;
        end
    end
    
    m  = 1;
    ID = find( strcmp( ar.model(m).data( dataSets(1) ).y, state ) );
    if ( isempty( ID ) )
        disp( 'Permissible states:' );
        fprintf( '%s\n', ar.model(m).data(dataSets(1)).y{:} );
        error( 'Could not find state %s', ID );
    end
    
    xExp     = log10( getConditionValues( ar, dataSets, variable, tp, [] ) );
    [~,i]    = sort( xExp );
    
    xSim     = log10( getConditionValues( ar, dataSets, variable, tp, 1 ) );
    [~,iSim] = sort( xSim );    
    
    yExp     = getFieldValues( ar, dataSets, 'tExp',  'yExp',           ID, tp, [] );
    yStd     = getFieldValues( ar, dataSets, 'tExp',  'ystdExpSimu',    ID, tp, [] );
    ySim     = getFieldValues( ar, dataSets, 'tFine', 'yFineSimu',      ID, tp, selectOne );
    ySimStd  = getFieldValues( ar, dataSets, 'tFine', 'ystdFineSimu',   ID, tp, selectOne );
    
    % For dose responses, the location of infinity needs to be placed manually
    if (nargin > 6)
        xExp( xExp == -inf ) = infinityLocation;
    end
    if (nargin > 6)
        xSim( xSim == -inf ) = infinityLocation;
    end    
    
    switch( errorBarMode )
        case 1
            h1          = errorbar( xExp(i), yExp(i), yStd(i), 'o' ); hold on;
            p1          = plot( xSim(iSim), ySim(iSim), '-' );
        case 2
            p1          = plot( xSim(iSim), ySim(iSim), '-' ); hold on;
            b1          = plot( xSim(iSim), ySim(iSim) + ySimStd(iSim) );
            b2          = plot( xSim(iSim), ySim(iSim) - ySimStd(iSim) );
            h1          = plot( xExp(i),    yExp(i), 'o' );
            set( b1, optsNoMarker{:}, 'LineWidth', .5 );
            set( b2, optsNoMarker{:}, 'LineWidth', .5 );
        case 3
            pat         = patch([xSim(iSim), fliplr(xSim(iSim))], [ySim(iSim) + ySimStd(iSim), fliplr(ySim(iSim) - ySimStd(iSim))], min(color + 0.4, [1,1,1]) ); hold on;
            p1          = plot( xSim(iSim), ySim(iSim), '-' );
            h1          = plot( xExp(i),    yExp(i), 'o' );
            set(pat, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.0, 'EdgeColor', [1,1,1]);
    end

    set( h1, opts{:}, 'LineWidth', 0.5 );
    set( p1, opts{:}, 'LineWidth', 1.0 );
    set( p1, 'Marker', 'none' );
    
    xmi = nanmin( xExp );
    xma = nanmax( xExp );
    ymi = nanmin( yExp - yStd );
    yma = nanmax( yExp + yStd );
    
    xL = xmi - abs((xma-xmi))*0.1;
    xU = xma + abs((xma-xmi))*0.1;
    yL = ymi - abs((yma-ymi))*0.1;
    yU = yma + abs((yma-ymi))*0.1;
    xlim( [xL, xU] );
    ylim( [yL, yU] );
end

% Get parameter values for a list of data sets
function cVals = getConditionValues( ar, datasets, conditionPar, tp, ID )
    cVals = [];
    
    for a = 1 : length( datasets )
        % User assumes this data set *should* only have one time point, so take that one
        if isempty( tp )
            if ( numel( unique(ar.model.data(datasets(a)).tExp ) ) > 1 )
                fprintf( '%f\n', ar.model.data(datasets(a)).tExp );
                error( 'Time point is ambiguous, please select a time point' );
            else
                tp = ar.model.data(datasets(a)).tExp(1);
            end
        end
        
        rep = numel(find(ar.model.data(datasets(a)).tExp == tp));

        if ( ~isempty(ID) )
            rep = 1;
        end
        
        for b = 1 : length( ar.model.data(datasets(a)).condition )
            if ( strcmp( ar.model.data(datasets(a)).condition(b).parameter, conditionPar ) )
                cVals = [ cVals repmat( str2num( ar.model.data(datasets(a)).condition(b).value ), 1, rep ) ]; %#ok
            end
        end
    end
end

function cVals = getFieldValues( ar, datasets, tField, field, input, tp, ID )

    cVals = [];
    
    pickByValue     = 0;
    if ( ( nargin < 4 ) || ( isempty( input ) ) )
        pickByValue     = 1;
        IDs             = 1;
    else
        if isnumeric( input )
            pickByValue = 1;
            IDs         = input;
        end
    end
    
    for a = 1 : length( datasets )
        if ( ~pickByValue )
            IDs = find( strcmp( ar.model.data(datasets(a)).y, input ) );
        end
        
        % User assumes this data set *should* only have one time point, so take that one
        if isempty( tp )
            if ( numel( unique(ar.model.data(datasets(a)).tExp ) ) > 1 )
                fprintf( '%f\n', ar.model.data(datasets(a)).tExp );
                error( 'Time point is ambiguous, please select a time point' );
            else
                tmp = ar.model.data(datasets(a)).tExp;
                tp = tmp(1);
            end
        end
        
        tmp = ar.model.data(datasets(a)).(tField);
        tID = find(tmp == tp);
        
        if ( ~isempty(ID) )
            tID = tID(ID);
        end
        
        if ( isempty( tID ) )
            fprintf( '%f\n', ar.model.data(datasets(a)).(tField) );
            error( 'Time point %g in data set %g does not exist, please select a time point', tp, datasets(a) );
        end

        cVals = [ cVals ; ar.model.data(datasets(a)).(field)(tID,IDs) ]; %#ok
    end
    cVals = cVals.';
end

