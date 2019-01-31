% arFprintf([output_level], [varargin] )
% 
% arFprintf wraps around fprintf. It uses global variable arOutputLevel to
% throw an output only, if output_level <= global arOutputLevel 
% 
% output_level <= arOutputLevel, then output is done 
%       0     - All printed output except errors suppressed
%       1     - Minimal output (only errors and requested info (arPrint etc))
%       2     - Regular output level (updates)
%       3     - Verbose output (everything)
% varargin    is passed to fprintf and contains format and output arguments
% 
% If arOutputLevel is not set, it is set to arOutputLevel = 2.
% A global variable is used to speed up the code. This function might be
% called very frequently.
% 
% Examples:
% >> global arOutputLevel
% >> arOutputLevel = 2;
% >> arFprintf(1,'Parameter %s = %d\n','kon',pi)
% Parameter kon = 3.141593e+00
% >> arFprintf(2,'Parameter %s = %.2f\n','kon',pi)
% Parameter kon = 3.14
% >> arFprintf(3,'Parameter %s = %.2f\n','kon',pi) % no output
% 



function arFprintf( output_level, varargin )
global arOutputLevel
if isempty( arOutputLevel )
    arOutputLevel = 2;
end

if ( output_level <= arOutputLevel )
    fprintf( varargin{:} );
end
