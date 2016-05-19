% Initialize and clear workspace of framework
%
% Data-2-Dynamics Software
% Website: https://bitbucket.org/d2d-development/d2d-software/wiki/Home
% Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de
% Copyright 2013 D2D Development Team. All rights reserved.

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
arFprintf(1, 'Copyright 2015 D2D Development Team. All rights reserved.\n\n');

ar.checksum = [];

if(~ispc)
    ar_path = strrep(which('arInit.m'),'/arFramework3/arInit.m','');
    [~,cmdout] = system(['hg summary -R ',ar_path]);
    leer = find(isspace(cmdout)==1);
    cmdout = cmdout(1:leer(2));
    ar.info.revision = cmdout;
    cmdout = [];
    leer = [];
else
    ar_path = strrep(which('arInit.m'),'\arFramework3\arInit.m','');
    [~,cmdout] = system(['hg summary -R ',ar_path]);
    leer = find(isspace(cmdout)==1);
    cmdout = cmdout(1:leer(2));
    ar.info.revision = cmdout;
    cmdout = [];
    leer = [];
end

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

clear j

