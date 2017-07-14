function [def_version_code, c_version_code] = arGetVersion

def_version_code = 3;
c_version_code = 'code_170714';

% Please note that if you add configuration flags that pertain to the C
% solver, please also add your configuration flag to arCheckCache, to make
% sure that the solver simulates the sensitivities when your settings are
% changed.
