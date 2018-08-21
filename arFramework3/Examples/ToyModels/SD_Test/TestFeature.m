function TestFeature()

global ar;

fprintf( 'TEST FOR STANDARD DEV PARAMETERS\n' );

arInit;
fprintf( 2, 'Parsing model... ' );
arLoadModel('SD_model');
arLoadData('ABC_data_BCobs');
arCompileAll(true);
fprintf( 'PASSED\n' );

fprintf( 2, 'Checking pError definition... ' );
arPrint   
if ( isequal(ar.qError,[0 1 0 0 0 1 1]) )
    fprintf('PASSED\n');
else
    error( 'WRONG ERROR ASSIGNMENT' );
end