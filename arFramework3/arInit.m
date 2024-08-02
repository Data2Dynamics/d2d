% arInit
% 
% Initialize and clear workspace of a new D2D project
% 
% Has to be called before defining a model and data in a script
%
% Example:
%   arInit;
%   arLoadModel('model1');
%   arLoadData('dataset1');
%   arCompileAll;
%   arFit;
%
% See also arLoadModel, arLoadData, arReset

global ar
arInitMain;


function arInitMain()

ar_path = fileparts(which('arInit.m'));
if(exist('arCheck','file') == 0)
    addpath([ar_path '/Subfunctions'])
end
% arCheck also adds paths of D2D subroutines. They have to be manually
% added there if a new folder is created in arFramework3.
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
arFprintf(1, 'Website & bug report: http://www.data2dynamics.org\n');
arFprintf(1, 'Copyright 2016 D2D Development Team. All rights reserved.\n\n');

ar.info.gitCommitHash = arGetGitCommitHash;

arCleanMemory;

ar.checksum = [];
ar.info.ar_path = ar_path;

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

% check toolbox requirements
verstr = ver;
req_toolboxes = {'Symbolic Math Toolbox', 'Optimization Toolbox'};
for i = 1:length(req_toolboxes)
    if ~ismember(req_toolboxes{i}, {verstr.Name})
        error('Missing toolbox: %s is required for using Data2Dynamics', req_toolboxes{i})
    end
end

ar.setup = struct;  % the setup commands are stored, here.
ar.setup.commands = cell(0);
ar.setup.arguments = cell(0);
ar.setup.commands{end+1} = mfilename; % this file name
ar.setup.arguments{end+1} = cell(0); % arInit has no arguments
ar.setup.modelfiles = {''};% model files to be read 
ar.setup.datafiles = {''}; % data files to be read 

% Get the setup file from the call stack
callStack = dbstack;
if length(callStack) > 2
    % arInit was called from another function/script -> probably setup file 
    ar.setup.setupfile = {callStack(3).file};
elseif length(callStack) == 2
    % arInit was called from the command line
    ar.setup.setupfile = {''};
end


ar = orderfields(ar);
ar.info = orderfields(ar.info);
ar.config = orderfields(ar.config);
ar.ppl = orderfields(ar.ppl);

end
