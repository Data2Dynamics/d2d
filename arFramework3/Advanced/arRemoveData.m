function arRemoveData(m,d,type,delayLink)
    
% arRemoveData(m,d,type,delayLink)
%
% Find and remove data added with arAddToData.
%
% m:         - Model index
% d:         - Data index
% type:      - Identifier which data points should be removed. This was set
%                 when calling arAddToData. [default:all]
% delayLink: - If set to 1, the model will not be relinked and arLink
%                 will have to be invoked.
%
% Note: Does this make arClearInterpolatedData obsolete?
%
% See also: arAddToData

    global ar;

    if ( nargin < 1 )
        help arRemoveData
        return
    end
    if ( ~isnumeric(m) || ( m > numel(ar.model) ) )
        error( 'Please specify a valid model index' );
    end
    if ( ( nargin < 2 ) || ~isnumeric(d) || d > numel( ar.model(m).data ) )
        error( 'Please specify a valid data index' );
    end
    if ( nargin < 3 )
        disp( 'Remove all added data points' );
        rm_id = find(ar.model(m).data(d).addedData);
    else
        rm_id = find(ar.model(m).data(d).addedData == type); 
    end
    if ( nargin < 4 )
        delayLink = 0;
    end
    
    ar.model(m).data(d).addedData(rm_id) = [];
    ar.model(m).data(d).tExp(rm_id) = [];
    ar.model(m).data(d).yExp(rm_id,:) = [];
    ar.model(m).data(d).yExpStd(rm_id,:) = [];    
    
    if ( delayLink ~= 1 )
        arLink;
    end

end

