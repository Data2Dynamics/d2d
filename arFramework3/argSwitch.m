% Function that puts arguments into an option struct
%
% Input arguments:
%   switches        -   Name of the arguments
%   extraArgs       -   Is the flag followed by an additional field with info?
%   description     -   Description of the flag (for output purposes)
%   verbose         -   Print output / Flags
%   varargin        -   varargin array which has the extra args

function [opts] = argSwitch( switches, extraArgs, description, verbose, varargin )

    for a = 1 : length(switches)
        opts.(lower(switches{a})) = 0;
        if ( extraArgs(a) == 0 )
            opts.([lower(switches{a}) '_args']) = {};
        end
    end
    
    a = 1;
    if ~isempty( varargin{1} )
        while( a <= length( varargin{1} ) )
            if (~ischar(varargin{1}{a}))
                error( 'Provided cell argument without indicating which switch it is for.' );
            end
            if ( max( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                error( sprintf( 'Invalid switch argument was provided "%s". Valid arguments are %s\n', varargin{1}{a}, sprintf( '''%s'' ', switches{:} ) ) ); %#ok
            end
            opts.(lower(varargin{1}{a})) = 1;
            if ( extraArgs( strcmpi( switches, varargin{1}{a} ) ) == 1 )
                try
                    opts.([lower(varargin{1}{a}) '_args']) = varargin{1}{a+1};
                    a = a + 1;
                catch
                    varargin{1}{a}
                    error( 'Did not provide arguments for flag %s', varargin{1}{a} );
                end
            end
            a = a + 1;
        end
    end
    if ( verbose )
        for a = 1 : length( switches )
            if ( extraArgs(a) == 0 )
                if ( ~isempty( description{a}{2-(opts.(lower(switches{a}))>0)} ) )
                    fprintf( '%s\n', description{a}{2-(opts.(lower(switches{a}))>0)} );
                end
            else
                if ( ~isempty( description{a}{1+(opts.(lower(switches{a}))>0)} ) )
                    fprintf( '%s\n', description{a}{1+(opts.(lower(switches{a}))>0)} );
                end
            end
        end
    end
    
end