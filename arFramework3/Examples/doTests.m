% Please verify any changes made to arCompile, arCompileAll, arLoadModel,
% arLoadData, arSimu, arSimuCalc, arFit by running these integration tests.

function doTests
    fprintf(2, 'This collection of tests checks whether specific functions\n');
    fprintf(2, 'in the D2D functions are working correctly. Please run this\n');
    fprintf(2, 'routine before pushing any changes to internal functions in D2D\n' );
    fprintf(2, 'to reduce the risk of pushing code that breaks existing\nfunctionality.\n\n' );

    try
        doTest('Advanced Events');
        doTest('Volume_Fitting_Test');
        doTest('MonotoneSpline');
        doTest('Stoichiometry');
    
        fprintf('Testing complete!\n');
    catch
        fprintf('Testing failed: please resolve errors and try again.\n');
        
    end
end

function doTest( dir )
    cd(dir);
    
    % Suppress output
    global arOutputLevel;
    arOutputLevel = 0;
    try
    TestFeature; cd('..');
    sprintf('\n');
    
    catch
        try
            fprintf(2, 'UNIT TEST FAILED: RESTARTING WITH VERBOSE OUTPUT\n\n' );
            arOutputLevel = 2;
            TestFeature;
        catch ME
            fprintf(getReport(ME));
            arOutputLevel = 0;
            cd('..');
            throw(MException('Testing:failedTest', 'Failed to pass a test. Please resolve the problem before pushing.'));
        end
    end
end