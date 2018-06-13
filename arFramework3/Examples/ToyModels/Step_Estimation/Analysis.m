% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. Since
% we are not optimizing over the location parameter, this is ok though.
%
% See estimateLocaton.m for an example where the location is also unknown.

setupStepInput

%%
arSimu(true,true,true); arChi2(true);
arPlotY; title('Pre-optimization (hit any key) ');  
pause;
arFit;
arSimu(true,true,true); arChi2(true);
arPlotY;


