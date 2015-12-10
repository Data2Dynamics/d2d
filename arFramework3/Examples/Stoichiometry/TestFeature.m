fprintf( 'INTEGRATION TEST FOR NON-UNITY STOICHIOMETRY\n' );

arInit(1);

fprintf( 2, 'Parsing model with non-unity stoichiometry ...\n' );
try
    arLoadModel('test2');
catch
    error('Test failed');
end

fprintf( 2, 'Compiling and simulating model with non-unity stoichiometry ...\n' );
try
    arCompileAll;
    arSimu(true, true, true);
catch
    error('Test failed');
end

Xi = ar.model.condition.xFineSimu(1,:);
X = ar.model.condition.xFineSimu(end,:);
Xall = ar.model.condition.xFineSimu;

fprintf( 2, 'Testing stoichiometry ratios for 3 X -> Y ... ' );
try
    if ( X( findState( ar, 'sB' ) ) / Xi( findState( ar, 'sA' ) ) < 0.33 ) || ...
        ( X( findState( ar, 'sD' ) ) / Xi( findState( ar, 'sC' ) ) > 0.34 )
        error('Test failed');
    else
        fprintf('PASSED\n');
    end
catch
    error('Test failed');
end

fprintf( 2, 'Testing stoichiometry ratio for X -> 3 Y ... ' );
try
    if ( X( findState( ar, 'sF' ) ) / Xi( findState( ar, 'sE' ) ) < 2.999 ) || ...
        ( X( findState( ar, 'sH' ) ) / Xi( findState( ar, 'sG' ) ) > 3.001 )
        error('Test failed');
    else
        fprintf('PASSED\n');
    end
catch
    error('Test failed');
end

fprintf( 2, 'Testing equality MASSACTION and CUSTOM ... ' );
try
    if ( norm( X( :, findState( ar, 'sA' ) ) - X( :, findState( ar, 'sC' ) ), 2 ) > 1e-6  ) || ...
        ( norm( X( :, findState( ar, 'sB' ) ) - X( :, findState( ar, 'sD' ) ), 2 ) > 1e-6  ) || ...
        ( norm( X( :, findState( ar, 'sE' ) ) - X( :, findState( ar, 'sG' ) ), 2 ) > 1e-6  ) || ...
        ( norm( X( :, findState( ar, 'sF' ) ) - X( :, findState( ar, 'sH' ) ), 2 ) > 1e-6  )
        error('Test failed');
    else
        fprintf('PASSED\n');
    end
catch
    error('Test failed');
end


