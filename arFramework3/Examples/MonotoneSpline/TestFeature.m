fprintf( 'INTEGRATION TEST FOR SPLINES\n' );

arInit;
fprintf( 2, 'Parsing model with cubic spline ...\n' );
try
    arLoadModel('normal_cubic');
    arLoadData('test', 1, 'csv');
    arLoadData('normal_cubic', 1, 'csv');
    arCompileAll(true);
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end

try
    fprintf( 2, 'Simulating and validating model with cubic spline ...\n' );
    arDisableData('normal_cubic');
    arFit;    
    if ( norm(ar.model.data(arFindData('normal_cubic')).res) < 1e-3 )
        fprintf('PASSED');
    end
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end

fprintf( 2, 'Parsing model with positive cubic spline ...\n' );
try
    arInit;
    arLoadModel('positive_cubic');
    arLoadData('test', 1, 'csv');
    arLoadData('positive_cubic', 1, 'csv');
    arCompileAll(2);
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end

try
    fprintf( 2, 'Simulating and validating model with cubic spline ...\n' );
    arDisableData('positive_cubic');
    arFit;    
    if ( norm(ar.model.data(arFindData('positive_cubic')).res) < 1e-3 )
        fprintf('PASSED');
    end
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end


fprintf( 2, 'Parsing model with monotonic spline ...\n' );
try
    arInit;
    arLoadModel('monotone');
    arLoadData('test', 1, 'csv');
    arLoadData('monotone', 1, 'csv');
    arCompileAll(2);
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end

try
    fprintf( 2, 'Simulating and validating model with cubic spline ...\n' );
    arDisableData('monotone');
    arFit;    
    if ( norm(ar.model.data(arFindData('monotone')).res) < 1e-3 )
        fprintf('PASSED');
    end
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end
