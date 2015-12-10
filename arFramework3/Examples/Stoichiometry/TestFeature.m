% Comparing MASSACTION with CUSTOM
arInit;
arLoadModel('test2');
arCompileAll;
arPlot;

Xi = ar.model.condition.xFineSimu(1,:);
X = ar.model.condition.xFineSimu(end,:);
Xall = ar.model.condition.xFineSimu;

fprintf( '> Integration test for non-unity stoichiometry\n' );

fprintf( 'Testing stoichiometry ratios for 3 X -> Y ... ' );
if ( X( findState( ar, 'sB' ) ) / Xi( findState( ar, 'sA' ) ) < 0.33 ) || ...
    ( X( findState( ar, 'sD' ) ) / Xi( findState( ar, 'sC' ) ) > 0.34 )
    error('Test failed');
else
    fprintf('PASSED\n');
end

fprintf( 'Testing stoichiometry ratio for X -> 3 Y ... ' );
if ( X( findState( ar, 'sF' ) ) / Xi( findState( ar, 'sE' ) ) < 2.999 ) || ...
    ( X( findState( ar, 'sH' ) ) / Xi( findState( ar, 'sG' ) ) > 3.001 )
    error('Test failed');
else
    fprintf('PASSED\n');
end

fprintf( 'Testing equality MASSACTION and CUSTOM ... ' );
if ( norm( X( :, findState( ar, 'sA' ) ) - X( :, findState( ar, 'sC' ) ), 2 ) > 1e-6  ) || ...
    ( norm( X( :, findState( ar, 'sB' ) ) - X( :, findState( ar, 'sD' ) ), 2 ) > 1e-6  ) || ...
    ( norm( X( :, findState( ar, 'sE' ) ) - X( :, findState( ar, 'sG' ) ), 2 ) > 1e-6  ) || ...
    ( norm( X( :, findState( ar, 'sF' ) ) - X( :, findState( ar, 'sH' ) ), 2 ) > 1e-6  )
    error('Test failed');
else
    fprintf('PASSED\n');
end


