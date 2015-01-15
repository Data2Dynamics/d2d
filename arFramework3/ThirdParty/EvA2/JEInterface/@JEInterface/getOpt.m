function val = getOpt(int, optName)
% Set a single optimset value within the JI instance.
% Arguments: 
%           int: the JEInterface instance
%           optName: name of the option to change, e.g. 'MaxFunEvals'
%           optVal: new value

val = int.opts.(optName);