
fprintf( 2, 'REGRESSION TEST FOR DOSE RESPONSE PREDICTOR FITTING\n' );

fprintf( 2, 'Loading model for dose response fitting ... ' );

% Load models & data
arInit;
arLoadModel('epo_binding');
arLoadData('Epo_binding_rep','epo_binding');
arCompileAll(true);

% parameters with prior information
arSetPars('init_Epo', log10(2100), 1, 1, -5, 4, 1, log10(2100), 0.02);

fprintf( 2, 'PASSED\n' );

arFit;
fprintf( 2, 'Testing final chi2fit ... ' );
if ( ar.chi2fit > -63.8033 )
    error( 'FINAL ERROR TOO LARGE' );
else
    fprintf( 2, 'PASSED\n' );
end
