% olist = arFindPar([ar], identifier, [varargin])
%
% Finds parameter(s) in the ar structure and returns their indices as a sorted vector.
% 
% ar            ar struct. Default is global [ar] struct.
% identifier    string or cell array of parameter names to be printed
% varargin      one or multiple optional strings, you can use to turn on the 
%               following options:
% 'names'       returns parameter names instead of indices
% 'verbose'     prints the parameter names to command window
% 'dynamic'     select only the dynamic parameters
% 'initial'     select only parameters contributing to initial values of states
% 'exact'       searches for exact match instead of seraching for pattern
% 'preserve'    output is kept in order of the input identifiers instead of
%               sorting it by indices.
%
% olist         list of indices/names in a (un)sorted order.
% 
% Examples
%   Return all parameter IDs containing "degrad" in the name
%       arFindPar( ar, 'degrad' )
%   Return the parameter ID corresponding to the parameter named degrad
%       arFindPar( ar, 'degrad', 'exact' )
%   Return all parameter IDs whose name contains "degrad" or "pro"     
%       arFindPar( ar, {'degrad', 'pro'} )
%   Return names of the parameters whose name contains "degrad" or "pro"
%       arFindPar( ar, {'degrad', 'pro'}, 'names' )
%   Show the names of the parameters it is returning
%       arFindPar( ar, {'degrad', 'pro'}, 'verbose' )
%   Return dynamic parameters only and returns them by name
%       arFindPar( ar, {'degrad', 'pro'}, 'dynamic', 'names' )
%   Return initials only and returns them by name
%       arFindPar( ar, {'degrad', 'pro'}, 'initial', 'names' )
%
% Note: If you wish to print formatted parameter values or get a parameter
% value, use arGetPars or arPrint instead.
% 
% See also  arGetPars arPrint

function olist = arFindPar( varargin )

    global ar;
    if ( isstruct( varargin{1} ) )
        arC = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    else
        arC = ar;
    end
    
    if ( ~iscell( varargin{1} ) )
        string{1} = varargin{1};
    else
        string = varargin{1};
    end

    opts = argSwitch( {'names', 'verbose', 'dynamic', 'initial', 'exact', 'preserve'}, varargin{2:end} );

    list = arC.pLabel;
    
    olist = cell(1, numel(string));
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
                if (~opts.dynamic || arC.qDynamic(a)) && ((~opts.initial || arC.qInitial(a)))
                    olist{b}(end+1) = a;
                end            
            end
        end
    end
    
    olist = [ olist{:} ];
    if ( ~opts.preserve )
        olist = unique( olist );
    end
    
    if ( opts.verbose )
        fprintf('%s\n', arC.pLabel{olist});
    end
    
    if ( opts.names )
        olist = arC.pLabel(olist);
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
    