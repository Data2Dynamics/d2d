fprintf( 'INTEGRATION TEST FOR VOLUME FITTING\n' );

arInit;
fprintf( 2, 'Parsing and compiling model with volume estimation... ' );
arLoadModel('test');
arLoadData('test', [], 'csv');
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Fitting model... ' );
ar.qFit=ones(size(ar.qFit));
arFit;
fprintf( 'PASSED\n' );

fprintf( 2, 'Testing whether fitting of volume was successful... ' );
if ( norm(ar.model.data.res) < 1e-1 )
    fprintf( 'PASSED\n' );
else
    error( 'FINAL ERROR TOO LARGE' );
end

fprintf( 2, 'Testing whether setting volume is successful... ' );
ar.p = [-0.6159, -0.3499, -0.4854];
arSimu(false, false, true);
if ( norm(ar.model.data.res) < 1e-1 )
    fprintf( 'PASSED\n' );
else
    error( 'FINAL ERROR TOO LARGE' );
end

fprintf( 2, 'Parsing and compiling model with volume estimation in the observation function... ' );
arInit;
arLoadModel( 'volUsedInDataDef' );
arLoadData( 'volUsedInDataDef', [], 'csv' );
arCompileAll(true);
ar.qFit=ones(size(ar.qFit));
fprintf( 'PASSED\n' );

fprintf( 2, 'Fitting model... ' );
arSimu(false,true,true);
arFit;
fprintf( 'PASSED\n' );

fprintf( 2, 'Testing whether fitting of volume was successful... ' );
estVol = arGetPars('vol_cyt', 0);
if ( (estVol > 9.9999 ) && (estVol < 10.0001) )
    fprintf( 'PASSED\n' );
else
    error( 'VOLUME ESTIMATE TOO FAR OFF' );    
end
if ( ar.chi2fit < 0.001 )
    fprintf( 'PASSED\n' );
else
    error( 'CHI2FIT TOO FAR OFF' );        
end
