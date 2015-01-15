function int = setOpt(int, optName, optVal)
% Set a single option value within the JI instance.
% Arguments: 
%           int: the JEInterface instance
%           optName: name of the option to change, e.g. 'MaxFunEvals'
%           optVal: new value

% this only makes sure the option actally exists and is valid
makeOptions(int, optName, optVal);

opts = int.opts;
opts.(optName) = optVal;
int.opts = opts;