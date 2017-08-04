arInit;
fprintf( 'INTEGRATION TEST FOR THE PREPROCESSOR\n' );
fprintf( 2, 'Loading model... ' );
arLoadModel('PreProc');
fprintf( 'PASSED\n' );

fprintf( 2, 'Checking if correct states have been included... ' );
if ( sum( ismember( ar.model.x, {'stateA', 'stateD', 'stateF', 'stateH', 'stateI', 'stateJ', 'stateL', 'stateO', 'stateR'} ) ) == 9 ) 
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end
fprintf( 2, 'Checking if correct states have been excluded... ' );
if ( sum( ismember( ar.model.x, {'stateB', 'stateC', 'stateE', 'stateG', 'stateK', 'stateM', 'stateN', 'stateP', 'stateQ'} ) ) == 0 )
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end

fprintf( 2, 'Loading data... ' );
arLoadData('PreProcData');
fprintf( 'PASSED\n' );

fprintf( 2, 'Checking if correct observables have been included... ' );
if ( sum( ismember( ar.model.data.y, {'dataA', 'dataD', 'dataF', 'dataH', 'dataI', 'dataJ', 'dataL', 'dataO', 'dataR'} ) ) == 9 ) 
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end
fprintf( 2, 'Checking if correct observables have been excluded... ' );
if ( sum( ismember( ar.model.data.y, {'dataB', 'dataC', 'dataE', 'dataG', 'dataK', 'dataM', 'dataN', 'dataP', 'dataQ'} ) ) == 0 )
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end
