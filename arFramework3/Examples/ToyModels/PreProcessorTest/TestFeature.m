function TestFeature()

global ar;

arInit;
fprintf( 'INTEGRATION TEST FOR THE PREPROCESSOR\n' );
fprintf( 2, 'Loading model... ' );
arLoadModel('PreProc');
fprintf( 'PASSED\n' );

fprintf( 2, 'Checking if correct states have been included... ' );
if ( sum( ismember( ar.model.x, {'stateA', 'stateD', 'stateF', 'stateH', 'stateI', 'stateJ', 'stateL', 'stateO', 'stateR', 'stateS', 'stateT', 'stateW', 'NESTER2', 'NESTEE3', 'stateZ', 'stateAfterComment' } ) ) == 16 ) 
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end
fprintf( 2, 'Checking if correct states have been excluded... ' );
if ( sum( ismember( ar.model.x, {'stateB', 'stateC', 'stateE', 'stateG', 'stateK', 'stateM', 'stateN', 'stateP', 'stateQ', 'stateU', 'stateV', 'stateX', 'stateY'} ) ) == 0 )
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end

fprintf( 2, 'Checking if named defines work (includes)...' );
if ( sum( ismember( ar.model.x, {'potato', 'PARA2'} ) ) == 2 ) 
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end
fprintf( 2, 'Checking if named defines work (excludes)...' );
if ( sum( ismember( ar.model.x, {'PARA1'} ) ) == 0 ) 
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

fprintf( 2, 'Checking if define blocks meant for replacement are parsing correctly... ' );
if ( strcmp(ar.model(end).fu{1}, '( ( this*is*a*test ) + ( ( f*a+p*b ) / ( 1 + pp * ( f*a+p*b ) ) ) )') )
    fprintf( 'PASSED\n' );
else
    error( 'FAILED\n');
end