% Initialize and clear workspace of framework
%
% Data-2-Dynamics Software
% Website: http://www.data2dynamics.org
% Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de
% Copyright 2016 D2D Development Team. All rights reserved.
ar_path = fileparts(which('arInit.m'));
if(exist('arCheck','file') == 0)
    addpath([ar_path '/Subfunctions'])
end
if(~arCheck)
    return;
end

global ar

warning('off', 'symbolic:sym:sym:DeprecateExpressions')

ar = struct([]);
ar(1).stop = 0;
ar.fevals = 0; 

arInitUser;

ar.info.initTime = now;
[ar.info.def_version_code, ar.info.c_version_code] = arGetVersion;
arFprintf(1, 'Data 2 Dynamics Software\n');
arFprintf(1, '(arFramework3, def-version %i, c-version %s)\n', ...
    ar.info.def_version_code, ar.info.c_version_code);
arFprintf(1, 'Website: http://www.data2dynamics.org\n');
arFprintf(1, 'Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de\n');
arFprintf(1, 'Copyright 2016 D2D Development Team. All rights reserved.\n\n');

ar.checksum = [];
ar.info.ar_path = ar_path;
clear ar_path;

% check for updates on github
arCheckVersion;

% initialize fields
ar = arInitFields(ar);

% check licenses
if ( ~license('test', 'Symbolic_Toolbox') )
    warning( 'D2D requires a license for the MathWorks symbolic math toolbox. It is unlikely that D2D will work.' );
end
if ( ~license('test', 'Optimization_Toolbox') )
    warning( 'No license found for optimization toolbox. If fitting is required, obtain a license or switch optimization method (e.g. ar.config.optimizer=3).' );
end

ar = orderfields(ar);
ar.info = orderfields(ar.info);
ar.config = orderfields(ar.config);
ar.ppl = orderfields(ar.ppl);
