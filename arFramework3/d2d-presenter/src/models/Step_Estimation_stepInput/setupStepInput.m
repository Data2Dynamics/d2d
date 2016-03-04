% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. Since
% we are not optimizing over the location parameter, this is ok though.
%
% See estimateLocaton.m for an example where the location is also unknown.

arInit;
arLoadModel('stepInput');
arLoadData('stepInput',1,'csv',true);

%% compile
arCompileAll(true);

ar.config.useEvents = 1;

% The position cannot be fitted with a normal step function
arSetPars('position', 50, 0, 0, 0, 100);

arFindInputs;
arLink;

% Deliberately set the parameter values to wrong values
arSetPars('before', 0.1, 1, 0, 0, 100);
arSetPars('after', 5, 1, 0, 0, 100);
arSetPars('degrad', .1  , 1, 0, 0, 100);

arSimu(true,true,true); arChi2(true);
%arPlotY; title('Pre-optimization (hit any key) ');  
%pause;
arFit;
arSimu(true,true,true); arChi2(true);
%arPlotY;


