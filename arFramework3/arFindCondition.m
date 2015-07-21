% Finds conditions in the ar structure and returns their condition indices 
% as a vector. Also checks for the specific condition parameters.
%
% Usage:
%   arFindCondition( ar, (model no), name, conditions )
% or:
%   arFindData( (ar), (model no), 'state', statename / number, conditions )
%
% Example:
%   arFindCondition( ar, (model no), 'mydata' )
%       Returns all condition IDs whose name contains "mydata"
%   arFindCondition( ar, (model no), 'mydata', 'dose', '100', 'actd', '1' )
%       Returns all condition IDs whose name contains "mydata" and who
%       correspond to dose being set to 100 and actd being set to 1.
%   arFindCondition( ar, (model no), 'mydata', 'verbose' )
%       Returns all condition IDs whose name contains "mydata" and prints
%       them.
%
% Examples of the second usage mode:
%    arFindData( ar, (model no), 'state', 'SOCS3' )
%       Returns all condition IDs that have the state SOCS3 in their observation
%       function.
%    arFindData( ar, (model no), 'state', 'SOCS3', 'verbose' )
%       Provides some additional debug information
%    arFindData( ar, (model no), 'state', 'SOCS3', 'names', 'il6', 100 )
%       Returns the same by name, and also filters on the value of il6 in
%       the experiment.
%
% The argument ar is optional. If not specified, the global ar structure is
% used. 
%
% Returns: List of IDs that correspond to the query.

function [olist] = arFindCondition( varargin )

    global ar;
    if ( isstruct( varargin{1} ) )
        ar = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    end
    
	[dlist, ~, m] = arFindData( varargin{:} );
    
    olist = [];
    for a = 1 : length( dlist )
        olist = union( olist, ar.model(m).data( dlist(a) ).cLink );
    end
end

