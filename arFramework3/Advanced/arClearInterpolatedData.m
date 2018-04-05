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

    if nargin < 1
        m = 1;
    end
    if nargin < 2
        ds = 'all';
    end
    
    all = 0;
    if ~isnumeric( ds )
        if strcmp( ds, 'all' )
            ds = 1 : numel( ar.model(m).data );
            all = 1;
        end
    end
    
    for jd = 1 : numel( ds )
        if ( ~isfield( ar.model(m).data(ds(jd)), 'addedData' ) )
            if ~all
                warning( 'Data %d does not have interpolated data', ds(jd) );
            end
        else
            interps = sum(ar.model(m).data(ds(jd)).addedData == 1, 2) > 0;
            
            if sum( interps ) > 0
                ar.model(m).data(ds(jd)).yExp(interps, :)    = [];
                ar.model(m).data(ds(jd)).yExpStd(interps, :) = [];
                ar.model(m).data(ds(jd)).tExp(interps)       = [];
                ar.model(m).data(ds(jd)).addedData(interps)  = [];
                fprintf( 'Deleted %d points from data ID %d\n', sum(interps), ds(jd) );
            else
                if ~all
                    warning( 'Data %d does not have interpolated data', ds(jd) );
                end
            end
        end
    end
    
    % Link the model
    arLink;
end
