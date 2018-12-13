% [olist] = arFindCondition( varargin )
%
% Finds conditions in the ar structure and returns their condition indices 
% as a vector. Also checks for the specific condition parameters.
%
% Possible inputs:
%   ar - arStruct
%   m - model index
%   name       - model specifiers like ('state',statename) or 
%                ('input',cell of input strings)
%   conditions - conditions
%
% Outputs:
%   olist - List of IDs that correspond to the query
%
% Usage:
%   arFindCondition( ar, (model no), name, conditions )
% or:
%   arFindCondition( (ar), (model no), 'state', statename / number, conditions )
% or
%   arFindCondition( (ar), (model no), 'input', cell array of input strings )
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
%    arFindCondition( ar, (model no), 'state', 'SOCS3' )
%       Returns all condition IDs that have the state SOCS3 in their observation
%       function.
%    arFindCondition( ar, (model no), 'state', 'SOCS3', 'verbose' )
%       Provides some additional debug information
%    arFindCondition( ar, (model no), 'state', 'SOCS3', 'names', 'il6', 100 )
%       Returns the same by name, and also filters on the value of il6 in
%       the experiment.
%
% Examples of the third usage mode:
%    arFindCondition( ar, (model no), 'input', {'input_dcf * step1(0, 10, 5)'} )
%       Returns all data IDs that have an input with the function
%       'input_dcf * step(0, 10, 5)'. Note that the input function is not
%       parsed and must appear exactly as in the .def file.
%
% By default arFindCondition is permissive. It will return all conditions that
% match the criterion and include those that don't match or where
% insufficient information is available whether they match. To reverse this 
% behaviour, add the flag 'conservative'. This will reject any dataset
% where the flag is not set.
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

