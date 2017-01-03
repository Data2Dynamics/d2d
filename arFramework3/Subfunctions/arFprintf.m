% arFprintf wraps around fprintf. It uses a global variable to determine
% the output level. The default arOutputLevel is 2.
%
% 0     - All printed output except errors suppressed
% 1     - Minimal output (only errors and requested info (arPrint etc))
% 2     - Regular output level (updates)
% 3     - Verbose output (everything)

function arFprintf( output_level, varargin )

    global arOutputLevel; 
    if isempty( arOutputLevel )
        arOutputLevel = 2;
    end

    if ( output_level <= arOutputLevel )
        fprintf( varargin{:} );
    end