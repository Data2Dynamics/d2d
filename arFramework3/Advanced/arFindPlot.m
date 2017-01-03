% Finds plots in the ar structure and returns their plot indices 
% as a vector.
%
% Usage:
%   arFindPlot( (ar), (model no), name, (flags) )
%
% Examples of the first usage mode:
%   arFindPlot( ar, (model no), 'myplot' )
%       Returns all plot IDs whose name contains "myplot"
%   arFindPlot( ar, (model no), {'myplot', 'potato'} )
%       Returns all plot IDs whose name contains "myplot" or "potato"
%   arFindData( ar, (model no), 'myplot', 'verbose' )
%       Returns all plot IDs whose name contains "myplot" and prints
%       them.
%   arFindData( ar, (model no), 'myplot', 'names' )
%       Only return the full names of the found plots
%
% Parameters enclosed by brackets are optional.
%
% Returns: List of IDs that correspond to the query and a cell array of the
% plot names which match the search criterion.

function [olist names m] = arFindPlot( varargin )

    if nargin == 0
        help arFindPlot;
        return;
    end

    global ar;
    if ( isstruct( varargin{1} ) )
        ar = varargin{1};
        if ( length( varargin ) > 1 )
            varargin = varargin(2:end);
        else
            error( 'Insufficient parameters' );
        end
    end
    
    % Did we specify a model number? If not, assume 1
    m = 1;
    if ( isnumeric( varargin{1} ) )
        m = varargin{1};
        if ( m > length(ar.model) )
            error( 'Model %d does not exist in ar structure (only has %d models)', m, length(ar.model) );
        end
        if (length( varargin ) > 1)
            varargin = varargin(2:end);
        else
            error( 'Need to supply data name' );
        end
    end
        
    if ( ischar( varargin{1} ) )
        string{1} = varargin{1};
        varargin = varargin(2:end);
    elseif ( iscell( varargin{1} ) )
        string = varargin{1};
        varargin = varargin(2:end);
    else
        error( 'Please supply a data name' );
    end

    switches    = { 'verbose', 'names', 'exact' };
    extraArgs   = [         0,       0,       0 ];
    
    opts = argSwitch( switches, extraArgs, varargin );
    
    verbose = opts.verbose;
    returnNames = opts.names;

    olist    = [];
    for b = 1 : length( string )
        for a = 1 : length( ar.model(m).plot )
            if ~opts.exact
                if ~isempty( findstr(lower(ar.model(m).plot(a).name), lower(string{b}) ) )
                    olist = union( olist, a );
                end
            else
                if strcmp(ar.model(m).plot(a).name, string{b})
                    olist = union( olist, a );
                end                
            end
        end
    end

    for a = 1 : length( olist )
        names{a} = ar.model(m).plot(olist(a)).name; %#ok<AGROW>
    end    
    
    if (verbose)
        for a = 1 : length( olist )
            fprintf( '%s\n', names{a} );
        end
    end
    
    if ( returnNames )
        olist = names;
    end
    
function [opts] = argSwitch( switches, extraArgs, varargin )

    for a = 1 : length(switches)
        if ( extraArgs(a) == 0 )
            opts.(lower(switches{a})) = 0;
        else
            opts.(lower(switches{a})) = {};
        end
    end
    
    a = 1;
    if ~isempty( varargin{1} )
        while( a <= length( varargin{1} ) )
            if ( max( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                error( 'Invalid switch argument was provided %s', varargin{1}{a} );
            end
            if ( extraArgs( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                opts.(lower(varargin{1}{a})) = 1;
            else
                try
                    opts.(lower(varargin{1}{a})) = varargin{1}{a+1};
                    a = a + 1;
                catch
                    error( 'Did not provide arguments for flag %s', varargin{1}{a} );
                end
            end
            a = a + 1;
        end
    end
    