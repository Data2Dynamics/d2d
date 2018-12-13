% arRemoveSplitPlots()
%
% Removes all run time generated plots (i.e. generated with arSplitPlot) from the model
%
% See also arSplitPlot
function arRemoveSplitPlots()
    global ar;
    
    for m = 1 : numel( ar.model )
        jps = [];
        for jp = 1 : numel( ar.model(m).plot )
            if ( isfield( ar.model(m).plot(jp), 'origin' ) && ~isempty( ar.model(m).plot(jp).origin ) )
                jps = [ jps, jp ];
            end
        end
        if ( numel( jps ) > 0 )
            arRemovePlot( m, jps );
            fprintf( 'Removed %d plots from model %d\n', numel( jps ), m );
        else
            fprintf( 'No plots were removed from model %d\n', m );
        end
    end
end