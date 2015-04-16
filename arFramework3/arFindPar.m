% Finds parameter(s) in the ar structure and returns their indices 
% as a vector.
%
% Usage:
%   arFindPar( (ar), names, (returnNames), (verbose) )
%
% Example:
%   arFindPar( ar, 'degrad' )
%       Returns all parameter IDs containing "degrad" in the name
%   arFindPar( ar, {'degrad', 'pro'} )
%       Returns all parameter IDs whose name contains "degrad" or "pro"
%   arFindPar( ar, {'degrad', 'pro'}, 'names' )
%       Returns names of the parameters whose name contains "degrad" or "pro"
%   arFindPar( ar, {'degrad', 'pro'}, 'verbose' )
%       Shows the names of the parameters it is returning
%
% The argument ar is optional. If not specified, the global ar structure is
% used.
%
% Note: If you wish to print formatted parameter values or get a parameter
% value, use arGetPars or arPrint instead.
%

function olist = arFindPar( varargin )

    global ar;
    if ( isstruct( varargin{1} ) )
        ar = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    end
    
    if ( ~iscell( varargin{1} ) )
        string{1} = varargin{1};
    else
        string = varargin{1};
    end
    
    names = 0;
    if ( length(varargin) > 1 )
        if ( strcmp( varargin{2}, 'names' ) )
            names = 1;
            if ( length( varargin ) > 1 )
                varargin = varargin(2:end);
            end
        end
    end    
    
    showNames = 0;
    if ( length(varargin) > 1 )
        if ( strcmp( varargin{2}, 'verbose' ) )
            showNames = 1;
        end
    end

    list = ar.pLabel;
    
    olist    = [];
    for b = 1 : length( string )
        for a = 1 : length( list )
             if ~isempty( strfind(lower(list{a}), lower(string{b}) ) )
                 olist = union( olist, a );
            end
        end
    end
    if ( showNames )
        fprintf('%s\n', ar.pLabel{olist});
    end
    
    if (names)
        olist = ar.pLabel(olist);
    end
    