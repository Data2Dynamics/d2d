% Function that puts arguments into an option struct
%
% Input arguments:
%   switches        -   Name of the arguments
%   extraArgs       -   Is the flag followed by an additional field with info?
%   description     -   Description of the flag (for output purposes)
%   verbose         -   Print output / Flags
%   varargin        -   varargin array which has the extra args
%
% Optionally, the first argument can be "softmatching" which matches known
% symbols but doesn't throw errors when unknown ones are found. The
% remaining arguments are left in varargin and returned.

function [opts, outargin] = argSwitch( switches, extraArgs, description, verbose, varargin )

    for a = 1 : length(switches)
        opts.(lower(switches{a})) = 0;
        if ( extraArgs(a) == 0 )
            opts.([lower(switches{a}) '_args']) = {};
        end
    end
    
    a = 1; soft = 0; outargin = {};
    if ~isempty( varargin{1} )
        % Soft matching, only respond to known arguments, and ignore the
        % others, removing the known ones off the stack
        if strcmp( varargin{1}, 'softmatching' )
            soft = 1;
            varargin = varargin(2:end);
        end
        while( a <= length( varargin{1} ) )            
            if (~ischar(varargin{1}{a}))
                if ( soft == 0 )
                    error( 'Provided cell argument without indicating which switch it is for.' );
                else
                    outargin{end+1} = varargin{1}{a}; %#ok
                end
            end
            if ( max( strcmpi( switches, varargin{1}{a} ) ) == 0 )
                if ( soft == 0 )
                    error( sprintf( 'Invalid switch argument was provided "%s". Valid arguments are %s\n', varargin{1}{a}, sprintf( '''%s'' ', switches{:} ) ) ); %#ok
                else
                    outargin{end+1} = varargin{1}{a};
                end
            else
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