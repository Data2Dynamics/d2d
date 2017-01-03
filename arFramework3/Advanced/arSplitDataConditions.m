% arSplitDataConditions( (ar), model, dataIDs, verbose )
%
% In the D2D system, conditions are typically used by multiple data files.
% When using the event system, this can be a problem; since one may want to
% set events for one data file, but not another (even though the
% corresponding condition is the same). This is where arSplitData comes in.
% When the condition associated with a dataset is shared with other datasets,
% arSplitData makes a copy of the conditions used in the dataset. This
% allows the user to set separate events for these newly created
% conditions.
%
% Usage:
%   arSplitDataConditions( model, dataIDs, verbose )
%
% Input arguments:
%   model       model number (typically 1)
%   dataIDs     data IDs (hint: you can find these using arFindData)
%   verbose     show what changes are made to the struct (1 is debug output on)  
%
% Note:
%   You should invoke this command *before* adding events. Otherwise, 
%   events will be copied along with the new conditions which can result 
%   in undefined behaviour.

function arSplitDataConditions( varargin )

    if ( nargin < 2 )
        error( 'Function needs at least two arguments.' );
    end
       
    global ar;
    
    % If we were called from a high level event function, since the higher
    % level function will already be in the command log and this one doesn't 
    % need to be added.
    s = dbstack(1);
    if ( (length(s)==0) || ( (~strcmp( s(1).file, 'arSteadyState.m' )) && (~strcmp( s(1).file, 'arFindInputs.m' )) ) )
        logCall( 'arSplitDataConditions', varargin{:} );
    end
    
    m = varargin{1};
    ds = varargin{2};
    if (nargin > 2)
        verbose = varargin{3};
    else
        verbose = 0;
    end
    
    for a = 1 : length( ds )
        d = ds(a);
        
        % Find the condition ID belonging to the dataset we wish to unique-ify
        c = ar.model(m).data(d).cLink;
        
        % Find which datasets this condition provides simulations for
        ds2 = ar.model(m).condition(c).dLink;
        
        if length(c) > 1
            error('Fatal error: Unexpected condition link length in ar.model(m).data(d).cLink');
        end
        
        % More than just the dataset we started with
        % ==> Needs to be unique-ified
        if ( length(ds2) > 1 )
            newID = length(ar.model(m).condition) + 1;
            % Add the new condition
            ar.model(m).condition(newID)        = ar.model(m).condition(c);
            ar.model(m).condition(newID).src    = c;
            ar.model(m).condition(newID).dLink  = d;
            
            % Remove the reference in the old data condition
            rem = find( ar.model(m).condition(c).dLink == d );
            ar.model(m).condition(c).dLink(rem) = [];
            
            % Add the reference to the newly created condition
            ar.model(m).data(d).cLink = newID;
            
            if ( verbose )
                fprintf( '[ C] Making new condition with condition ID %d for data %s with data ID %d\n', newID, ar.model(m).data(d).name, d );
            end
        else
            if ( verbose )
                % This one is already unique
                fprintf( '[NC] Condition corresponding to data %s with data ID %d is already unique\n', ar.model(m).data(d).name, d );
            end
        end
    end
    
end

% Event logging
function logCall( fn, varargin )
    global ar;
    
    if ~isfield(ar, 'eventLog')
        ar.eventLog = {};
    end
    call = [fn '('];
    if length( varargin ) > 0
        call = sprintf('%s%s', call, mat2str(varargin{1}) );
    end
    for a = 2 : length( varargin )
        call = sprintf('%s, %s', call, mat2str(varargin{a}) );
    end
    call = [call ')'];
    ar.eventLog{length(ar.eventLog)+1} = call;
end