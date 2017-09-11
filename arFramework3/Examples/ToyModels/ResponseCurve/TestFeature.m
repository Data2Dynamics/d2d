function TestFeature()

fprintf( 2, 'INTEGRATION TEST FOR BI-LINEAR LOOKUP TABLE (LUT)\n' );
fprintf( 2, 'Loading model for LUT test... ' );

arInit;
arLoadModel('responseCurve');
arLoadData('responseData', 1, 'csv');
arCompileAll(true);
ar.config.rtol=1e-9; ar.config.atol=1e-9;

% The dynamics of the system are set so fast that the state will follow the
% response curve exactly (equilibrate to it rapidly)
ar.p(1) = 10;
ar.ub(1) = 10;

% Make sure we can get the right points for testing the surface
tDesired = 0: 1/4 : 1;
ar.model(1).data(1).tExtra = tDesired;
arLink;

fprintf( 'PASSED\n' );
fprintf( 2, 'Testing LUT outcomes... ' );

for c = 0 : 8
    arSetPars( 'parameter', (c/8), 1, 0 ); arSimu(false, true, true);
    p(c+1, :) = ar.model(1).condition(1).xFineSimu( ismember( ar.model(1).condition(1).tFine, tDesired ), 1 ); %#ok
end
trueValues = [1,2,3,4,5; 1.5,2.5,3.5,4.5,5.5; 2,3,4,5,6; 2.5,3.5,4.5,5.5,6.5; 3,4,5,6,7; 3.5,4.5,5.5,6.5,7.5; 4,5,6,7,8; 4.5,5.5,6.5,7.5,8.5; 5,6,7,8,9];

if ( sum(sum((trueValues-p).^2)) < 1e-7 )
    fprintf('PASSED\n');
else
    error( 'FAILED, WRONG INTERPOLATION IN LUT' );
end

fprintf( 2, 'Testing LUT sensitivities... ' );
arSetPars( 'parameter', .56, 1, 0 ); arSimu(false, true, true);
arCalcMerit;
arFiniteDifferences(1e-6);
dif = ar.sresFD - ar.sres;
diffmerge = zeros( size( dif ) );
for i=1:size(dif,1), for j=1:size(dif,2), diffmerge(i,j)=min([dif(i,j)/ar.sres(i,j), dif(i,j)]); end, end;
if ( max(max(diffmerge)) < 1e-5 )
    fprintf('PASSED\n');
else
    error( 'FAILED, SENSITIVITY ERROR TOO LARGE\n' );
end