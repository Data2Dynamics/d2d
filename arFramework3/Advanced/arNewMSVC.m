% arNewMSVC
%
% Adds a compilation flag which avoids redefinition of TIMESPEC during compilation
% 
% ** only for Microsoft Visual C++ **
%
% See wiki page: https://github.com/Data2Dynamics/d2d/wiki/I-get-the-error-%22Error-C2011-'timespec':-'struct'-type-redefinition%22-while-compiling


function arNewMSVC()

    arAddCompilationFlag( '-DHAVE_STRUCT_TIMESPEC' );
    
    