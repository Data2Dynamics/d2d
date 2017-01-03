% Finds parameter(s) in the ar structure and returns their indices 
% as a vector.
%
% Usage:
%   arFindPar( (ar), names, (returnNames), (verbose), (dynamic), (exact) )
%
% Example:
%   arFindPar( ar, 'degrad' )
%       Returns all parameter IDs containing "degrad" in the name
%   arFindPar( ar, 'degrad', 'exact' )
%       Returns the parameter ID corresponding to the parameter named degrad
%   arFindPar( ar, {'degrad', 'pro'} )
%       Returns all parameter IDs whose name contains "degrad" or "pro"
%   arFindPar( ar, {'degrad', 'pro'}, 'names' )
%       Returns names of the parameters whose name contains "degrad" or "pro"
%   arFindPar( ar, {'degrad', 'pro'}, 'verbose' )
%       Shows the names of the parameters it is returning
%   arFindPar( ar, {'degrad', 'pro'}, 'dynamic', 'names' )
%       Only returns dynamic parameters and returns them by name
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
        arOld = ar;
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

    opts = argSwitch( {'names', 'verbose', 'dynamic', 'exact'}, varargin{2:end} );

    list = ar.pLabel;
    
    olist    = [];
    for b = 1 : length( string )
        for a = 1 : length( list )
            if ( opts.exact )
                if strcmp(lower(list{a}), lower(string{b}) )
                    l = 1;
                else
                    l = [];
                end
            else
                l = strfind(lower(list{a}), lower(string{b}) );
            end
            if ~isempty( l )
                if (~opts.dynamic || ar.qDynamic(a))
                    olist = union( olist, a );
                end
            end
        end
    end
    if ( opts.verbose )
        fprintf('%s\n', ar.pLabel{olist});
    end
    
    if ( opts.names )
        olist = ar.pLabel(olist);
    end
    
    if ( exist( 'arOld' ) )
        ar = arOld;
    end
end

function [opts] = argSwitch( switches, varargin )

    for a = 1 : length(switches)
        opts.(switches{a}) = 0;
    end

    for a = 1 : length( varargin )
        if ( max( strcmp( varargin{a}, switches ) ) == 0 )
            error( 'Invalid switch argument was provided %s', varargin{a} );
        end
        opts.(varargin{a}) = 1;
    end    
end
    