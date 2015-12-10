fprintf( 'INTEGRATION TEST FOR VOLUME FITTING\n' );

arInit;
fprintf( 2, 'Parsing model with volume estimation ...\n' );
try
    arLoadModel('test');
    arLoadData('test', [], 'csv');
    arCompileAll(true);
catch ME
    fprintf(getReport(ME));
    error('FAILED');
end

fprintf( 2, 'Testing whether fitting of volume was successful ...' );
if ( norm(ar.model.data.res) < 1e-1 )
    disp( 'PASSED' );
else
    error( 'FAILED' );
end