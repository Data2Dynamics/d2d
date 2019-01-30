% [def_version_code, c_version_code] = arGetVersion
%
% Returns the version of the definition file and c version of the code.
% Note that whenever changes are made in the C file, this version number
% needs to be updated to trigger a recompile on other models.
%
%   def_version_code    Def format version
%   c_version_code      C version of the code
function [def_version_code, c_version_code] = arGetVersion

def_version_code = 3;
c_version_code = 'code_180926b';

% Please note that if you add configuration flags that pertain to the C
% solver, please also add your configuration flag to arCheckCache, to make
% sure that the solver simulates the sensitivities when your settings are
% changed.
