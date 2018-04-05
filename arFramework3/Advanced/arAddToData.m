% Function which adds the data points to the arStruct
%
% Usage:
%   m           - Model idx
%   d           - Data idx (find with arFindData)
%   obs         - Observable idx or name to add data for
%   tExp        - New time points to add
%   yExp        - y values associated with observable obs for the new timepoints
%   type        - Type identifier (integer which helps identify added data)
%                   1 is reserved for interpolated data from arInterpolateData
%   delayLink   - If set to 1, the model will not be relinked and arLink
%                 will have to be invoked.
function arAddToData( m, d, obs, tExpNew, yExpNew, type, delayLink )
    global ar;

    if ( nargin < 1 )
        help arAddToData
        return
    end
    if ( ~isnumeric(m) || ( m > numel(ar.model) ) )
        error( 'Please specify a valid model index' );
    end
    if ( ( nargin < 2 ) || ~isnumeric(d) || d > numel( ar.model(m).data ) )
        error( 'Please specify a valid data index' );
    end
    if ( nargin < 3 )
        error( 'Please specify an observable' );
    end
    if ( nargin < 4 )
        error( 'Please specify time points' );
    end
    if ( nargin < 5 )
        error( 'Please specify observable points' );
    end
    if ( nargin < 6 )
        error( 'Please specify an identifier for the new data' );
    end
    if ( type == 0 )
        error( 'Type == 0 is not allowed, since this would not allow the system to identify this added data later' );
    end
    if ( numel( tExpNew ) ~= numel( yExpNew ) )
        error( 'Dimensions of tExp and yExp must be the same' );
    end
    if ( nargin < 7 )
        delayLink = 0;
    end
    if ( isnumeric( obs ) )
        if ( obs > numel(ar.model(m).data(d).y) )
            error( 'Invalid observable index' );
        end
    else
        obs = find( strcmp( obs, ar.model(m).data(d).y ) );
        if ( isempty( obs ) )
            error( 'Invalid observable. Available ones for that dataset are: %s', sprintf('%s ', ar.model(m).data(d).y{:}) );
        end
    end
    
    nExps = numel( yExpNew );
    nans = nan( nExps, numel( ar.model(m).data(d).y ) );
    nData = nans;
    nData(:, obs) = yExpNew;    
    
    if ~isfield( ar.model(m).data(d), 'addedData' )
        ar.model(m).data(d).addedData = [];
    end
    
    ar.model(m).data(d).addedData( size( ar.model(m).data(d).yExp, 1 ) + 1 : size( ar.model(m).data(d).yExp, 1 ) + size( nData, 1 ), : ) = type;
    ar.model(m).data(d).yExp = [ ar.model(m).data(d).yExp; nData ];
    ar.model(m).data(d).tExp = [ ar.model(m).data(d).tExp; tExpNew.' ];
    ar.model(m).data(d).yExpStd = [ ar.model(m).data(d).yExpStd; nans ];
    
    if ( delayLink ~= 1 )
        arLink;
    end
end