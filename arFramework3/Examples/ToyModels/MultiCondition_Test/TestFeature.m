function TestFeature()

% This file tests optimization of models with complex step functions
% It is expected to throw a warning about the location parameter. Since
% we are not optimizing over the location parameter, this is ok though.
fprintf( 'INTEGRATION TEST FOR RANDOM PARAMETER EFFECTS\n' );

fprintf( 2, 'Loading model for random parameter effect test... ' );
arInit;
arLoadModel('stepInput');
arLoadData('stepInput',1,'csv',true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Compiling model for step location test... ' );
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Estimating ... ' );
arFit;
arFit;
arFit;
fprintf( 'PASSED\n' );

fprintf( 2, 'Determining whether fit is acceptable... ' );
if (ar.chi2 < 1e-4)
    fprintf( 'PASSED\n' );
else
    error( 'FINAL ERROR TOO LARGE' );
end