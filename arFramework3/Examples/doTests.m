% Please verify any changes made to arCompile, arCompileAll, arLoadModel,
% arLoadData, arSimu, arSimuCalc, arFit by running these integration tests.

function doTests
    try
        doTest('Advanced Events');
        doTest('Volume_Fitting_Test');
        doTest('MonotoneSpline');
        doTest('Stoichiometry');
    catch
    end
    fprintf( 'All tests passed!\n' );
end

function doTest( dir )
    cd(dir);
    
    % Suppress output
    global arOutputLevel;
    arOutputLevel = 0;
    try
    TestFeature;
    sprintf('\n');
    catch
        try
            fprintf(2, 'UNIT TEST FAILED: RESTARTING WITH VERBOSE OUTPUT\n\n' );
            arOutputLevel = 2;
            TestFeature;
        catch ME
            fprintf(getReport(ME));
            arOutputLevel = 0;
            rethrow('Failed to pass a test. Please resolve the problem before pushing.');
        end
    end
    cd('..');
end