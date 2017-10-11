function qLabels( lab, xlabeln, ylabeln, xt, yt, N, N2, logtrafo, infinityLocation )

    if (nargin < 6)
        N = 1;
    end
    if (nargin < 7)
        N2 = 1;
    end
    if ( nargin < 8)
        logtrafo = [0,0];
    end


    
    box off;
    
    set( gca, 'LineWidth', 0.5 );
        if ( ~isempty(xt) )
            xr = traf(xt, logtrafo(1));

            % There is a -inf in the list; replace it!
            if ( max( xr == -inf ) > 0 )
                xr( xr == -inf ) = min( xr(xr~=-inf) ) - 1;
            end

            % Set x limits
            d = 0.1*(max(xr)-min(xr));
            xlim([min(xr)-d, max(xr)+d]);            
            for a = 1 : length( xt )
                xtl{a} = xt(a); %#OK<AGROW>
            end            
            set( gca, 'XTickMode', 'manual' );
            set( gca, 'XTick', xr );
            set( gca, 'XTickLabel', xtl );
        end
        if ( ~isempty(yt) )
            for a = 1 : length( yt )
                ytl{a} = yt(a); %#OK<AGROW>
            end   
            set( gca, 'YTickMode', 'manual' );    
            set( gca, 'YTickLabel', ytl );
            set( gca, 'YTick', traf(yt, logtrafo(2)) );
        end
 
    set( gca, 'TickDir', 'out' );
    set( gca, 'TickLength', [0.025, 0.025] );
    %set( gca, 'Position', [0.1 + (N-1) * 0.3, 0.4 + (N2-1)*0.25, 0.2, 0.15] );
    %set( gca, 'Position', [0.05 + (N-1) * 0.31, 0.4 + (N2-1)*0.26, 0.2, 0.15] );
    global errorBarMode
    if ( errorBarMode < 3 )
        set( gca, 'Position', [0.1 + (N-1) * 0.3, 0.4 + (N2-1)*0.25, 0.2, 0.15] );
        %set( gca, 'Position', [0.05 + (N-1) * 0.31, 0.1 + (N2-1)*0.35, 0.2, 0.15] );
    else
        dalt = 0.05;
        set( gca, 'Position', [-0.2 + (N-1) * (0.33+dalt), 0 + (N2-1)*(0.35+dalt), 0.2+dalt, 0.15+dalt] );
    end
    title( lab );
    
    set(gcf, 'Clipping', 'off' );
    set(gcf, 'Renderer', 'painters' );
    set(gcf, 'RendererMode', 'manual' );
    set(gcf, 'PaperType', '<custom>' );
    set(gcf, 'PaperUnits', 'centimeters' );
    set(gcf, 'PaperSize', [20, 20] )
    set(gcf, 'Position', [0 0 20 20] );
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 20 20] );
    set(gcf, 'Resize', 'On' );
    
    ylabel(ylabeln);
    xlabel(xlabeln);
end

function v = traf( v, trafo )
    switch( trafo )
        case 1
            v = log10( v );
        case -1
            v = 10.^v;
        case 0
        otherwise
            error('Invalid transform');
    end
end
