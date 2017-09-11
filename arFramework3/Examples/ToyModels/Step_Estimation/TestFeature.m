function TestFeature()

% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. Since
% we are not optimizing over the location parameter, this is ok though.
fprintf( 'INTEGRATION TEST FOR STEP LOCATION ESTIMATION\n' );

fprintf( 2, 'Loading model for step location test... ' );
arInit;
arLoadModel('smoothStepInput');
arLoadData('stepInput',1,'csv',true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Compiling model for step location test... ' );
arCompileAll(true);
fprintf( 'PASSED\n' );

% In the smooth step setting it can be fitted without problems
arSetPars('position', 22, 1, 0, 0, 100);

% If we look at the data, we can see that the time over which the step is
% active is about 25 time units long. We add some events on this range to
% make sure that the integrator does not step over the step entirely.
arAddEvent(1, 'All', 0:25:100);

% The problem with step functions, is that their 'change' is only in a very
% limited range. In order to find the right time point, the stepping
% behaviour should be wide in a first fitting step (where the location is
% initially determined, hence we add a lower bound on the smoothness.
arSetPars('smoothness', 0, 1, 1, 0, 3);

% Deliberately set the parameter values to wrong values
arSetPars('after', 5, 1, 0, 0, 100);
arSetPars('degrad', .1, 1, 0, 0, 100);

fprintf( 2, 'Estimating phase I... ' );
arFit;
fprintf( 'PASSED\n' );

fprintf( 2, 'Estimating phase II... ' );
arSetPars('smoothness', 0, 1, 1, -5, 5);
arFit;
arFit;
arFit;
fprintf( 'PASSED\n' );

fprintf( 2, 'Determining whether fit is acceptable... ' );
if (norm(ar.model.data.yExp-ar.model.data.yExpSimu)<1e-3)
    fprintf( 'PASSED\n' );
else
    error( 'FINAL ERROR TOO LARGE' );
end