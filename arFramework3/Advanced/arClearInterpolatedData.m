% Remove datapoints which were added with arInterpolateData
%
% Usage:
%
%   function arClearInterpolatedData( m, ds, obs, tmin, tmax )
%
%    m        - Model index
%    ds       - Data indices, find using arFindData ('all' also works)
%

function arClearInterpolatedData( m, ds )  

    global ar;

    if nargin < 2
        help arClearInterpolatedData;
        error( 'Function needs at least 2 arguments.' );
    end
    
    if ~isnumeric( ds )
        if strcmp( ds, 'all' )
            ds = 1 : numel( ar.model(m).data );
            all = 1;
        end
    end
    
    for jd = 1 : numel( ds )
        if ( ~isfield( ar.model(m).data(ds(jd)), 'interpolatedData' ) )
            if ~all
                warning( 'Data %d does not have interpolated data' );
            end
        else
            interps = sum(ar.model(m).data(ds(jd)).interpolatedData, 2) > 0;
            
            if sum( interps ) > 0
                ar.model(m).data(ds(jd)).yExp(interps, :)    = [];
                ar.model(m).data(ds(jd)).yExpStd(interps, :) = [];
                ar.model(m).data(ds(jd)).tExp(interps)       = [];
                ar.model(m).data(ds(jd)).interpolatedData    = [];
                fprintf( 'Deleted %d points from data ID %d\n', sum(interps), ds(jd) );
            else
                if ~all
                    warning( 'Data %d does not have interpolated data' );
                end
            end
        end
    end
    
    % Link the model
    arLink;
end
