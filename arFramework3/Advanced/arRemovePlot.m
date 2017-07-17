% Remove a plot from a specific model
% 
% Usage:
%   arRemovePlot( m, p )
%
%   m       - Model ID
%   ps      - Plot IDs (0 indicates last plot)
%
% Hint: Determine plot IDs by arFindPlot

function arRemovePlot(m, ps)
    global ar;
    
    if ( m > numel( ar.model ) )
        error( 'Model does not exist' );
    end
    
    for jp = 1 : numel( ps )
        if ( ps(jp) < 1 )
            ps(jp) = numel( ar.model(m).plot ) - p(jp);
        end
        if ( ps(jp) > numel( ar.model(m).plot ) )
            error( 'Plot %d does not exist', p(jp) );
        end
    end
    
    ar.model(m).plot(ps) = [];
    ar.model(m).qPlotXs(ps) = [];
    ar.model(m).qPlotYs(ps) = [];
    ar.model(m).qPlotVs(ps) = [];
end