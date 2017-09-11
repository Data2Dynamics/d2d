function TestFeature()

global ar;

fprintf( 'INTEGRATION TEST FOR FLUX FITTING\n' );

arInit;
fprintf( 2, 'Parsing and compiling model with flux fitting... ' );

arInit;
arLoadModel('test');
arLoadData('test', [], 'csv');
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Fitting model... ' );
arFit;

estFlux = ar.model.condition.vExpSimu(end,ismember(ar.model.v, 'Fittedflux'));

fprintf( 2, 'Checking flux... ' );
if ( ( estFlux > 1.99999 ) && ( estFlux < 2.00001 ) )
    fprintf( 'PASSED\n' );
else
    error( 'FINAL DEVIATION TOO LARGE' );
end