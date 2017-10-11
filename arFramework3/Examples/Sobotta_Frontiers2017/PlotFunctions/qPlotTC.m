function h1 = qPlotTC( ar, data, color, ID )

    global errorBarMode;
    
    opts         = {'Color', color, 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'Marker', 'o', 'MarkerSize', 3};
    optsNoMarker = {'Color', min(color + 0.4, [1,1,1])};
    
    m = 1;
    try
        ID      = find( strcmp( ar.model(m).data( data ).y, ID ) );
    catch
        ar.model(m).data( data ).y
        ID                                          %#ok
        data                                        %#ok
        strcmp( ar.model(m).data( data ).y, ID )
    end
    if ( isempty( ID ) )
        disp( 'Permissible states:' );
        fprintf( '%s\n', ar.model(m).data(data).y{:} );
        error( 'Could not find state %s', ID );
    end
    
    tExp    = ar.model(m).data( data ).tExp;
    tFine   = ar.model(m).data( data ).tFine;
    yExp    = ar.model(m).data( data ).yExp(:, ID);
    yFine   = ar.model(m).data( data ).yFineSimu(:, ID);
    ystd    = ar.model(m).data( data ).ystdExpSimu(:, ID);
    yfinestd= ar.model(m).data( data ).ystdFineSimu(:, ID);
    
    switch( errorBarMode )    
        case 1
            h1          = errorbar( tExp, yExp, ystd, 'o' ); hold on;
            p1          = plot( tFine, yFine, '-' );
        case 2
            p1          = plot( tFine, yFine, '-' ); hold on;
            b1          = plot( tFine, yFine + yfinestd );
            b2          = plot( tFine, yFine - yfinestd );
            h1          = plot( tExp,  yExp,  'o' );
            set( b1, optsNoMarker{:}, 'LineWidth', .5 );
            set( b2, optsNoMarker{:}, 'LineWidth', .5 );
        case 3
            pat         = patch([tFine; flipud(tFine)], [yFine+yfinestd; flipud(yFine-yfinestd)], min(color + 0.4, [1,1,1]) ); hold on;
            p1          = plot( tFine, yFine, '-' );
            h1          = plot( tExp,  yExp,  'o' );
            set(pat, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.0, 'EdgeColor', [1,1,1]);
    end

    set( h1, opts{:}, 'LineWidth', 0.5 );
    set( p1, opts{:}, 'LineWidth', 1.0 );
    set( p1, 'Marker', 'none' );
    
    xmi = min( tFine );
    xma = max( tFine );
    ymi = (min( [yFine - yfinestd; yExp - ystd] ));
    yma = (max( [yFine + yfinestd; yExp + ystd] ));
    
    xL = xmi - abs((xma-xmi))*0.1;
    xU = xma + abs((xma-xmi))*0.1;
    yL = ymi - abs((yma-ymi))*0.1;
    yU = yma + abs((yma-ymi))*0.1;

    xlim( [xL, xU] );
    ylim( [yL, yU] );
end


