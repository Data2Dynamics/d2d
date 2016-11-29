% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. However,
% this is not a problem in the case of this function.
% 
% These functions interpolate between 0 and 1. The given parameters indicate 
% the start and end of the transition period. Smooth1 is first order continuous, 
% while smooth2 is second order continuous.

arInit;
arLoadModel('smoothies');
arLoadData('smoothies',1,'csv',true);

%% compile
arCompileAll(false);
ar.config.showFitting   = 1;

ar.config.optim.TolFun  = 1e-4;
ar.config.optim.TolX    = 1e-4;

% Put start and end at both extremes of the time course
ar.p( arFindPar( 't1' ) ) = log10( 0.001 );
ar.p( arFindPar( 't2' ) ) = log10( 100 );

arSimu(true,true,true); arChi2(true);
arPlotY; title('Pre-optimization (hit any key) ');  
pause;
title('Optimizing (step 1)');  
arFit;

title('Optimizing (step 2)');  
arFit;