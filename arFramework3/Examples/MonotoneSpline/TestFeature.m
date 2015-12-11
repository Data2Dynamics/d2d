fprintf( 'INTEGRATION TEST FOR SPLINES\n' );

arInit;
fprintf( 2, 'Parsing model with cubic spline... ' );
arLoadModel('normal_cubic');
arLoadData('test', 1, 'csv');
arLoadData('normal_cubic', 1, 'csv');
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Simulating and validating model with cubic spline... ' );
arDisableData('normal_cubic');
arFit;    
if ( norm(ar.model.data(arFindData('normal_cubic')).res) < 3e-3 )
    fprintf('PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end

fprintf( 2, 'Parsing model with positive cubic spline... ' );
arInit;
arLoadModel('positive_cubic');
arLoadData('test', 1, 'csv');
arLoadData('positive_cubic', 1, 'csv');
arCompileAll(2);
fprintf( 'PASSED\n' );
fprintf( 2, 'Simulating and validating model with positive cubic spline... ' );
arDisableData('positive_cubic');
arFit;    
if ( norm(ar.model.data(arFindData('positive_cubic')).res) < 3e-3 )
    fprintf('PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end


fprintf( 2, 'Parsing model with monotonic spline... ' );
arInit;
arLoadModel('monotone');
arLoadData('test', 1, 'csv');
arLoadData('monotone', 1, 'csv');
arCompileAll(2);
fprintf( 'PASSED\n' );

fprintf( 2, 'Simulating and validating model with cubic spline... ' );
arDisableData('monotone');
arFit;    
if ( norm(ar.model.data(arFindData('monotone')).res) < 3e-3 )
    fprintf('PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end
