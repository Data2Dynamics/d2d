arInit;
arLoadModel('responseCurve');
arLoadData('responseData', 1, 'csv');
arCompileAll(true);

% The dynamics of the system are set so fast that the state will follow the
% lookup table exactly (equilibrate to it rapidly)
ar.p(1) = 10;
ar.ub(1) = 10;

% This stuff needs pretty tight tolerances
ar.config.rtol=1e-9;
ar.config.atol=1e-9;

% Make sure we can get the right points for testing the surface
tDesired = 0: 1/4 : 1;
ar.model(1).data(1).tExtra = tDesired;
arLink;

for c = 0 : 8
    arSetPars( 'parameter', (c/8), 1, 0 ); arSimu(false, true, true);
    p(c+1, :) = ar.model(1).condition(1).xFineSimu( ismember( ar.model(1).condition(1).tFine, tDesired ), 1 ); %#ok
end
trueValues = [1,2,3,4,5; 1.5,2.5,3.5,4.5,5.5; 2,3,4,5,6; 2.5,3.5,4.5,5.5,6.5; 3,4,5,6,7; 3.5,4.5,5.5,6.5,7.5; 4,5,6,7,8; 4.5,5.5,6.5,7.5,8.5; 5,6,7,8,9];
imagesc( p );
title( 'ResponseCurve' );

figure; arPlot;
