% Please verify any changes made to arCompile, arCompileAll, arLoadModel,
% arLoadData, arSimu, arSimuCalc, arFit by running these integration tests.

function doTests
    global ar;

    try
        doTest('Advanced Events');
        doTest('Volume_Fitting_Test');
        doTest('MonotoneSpline');
        doTest('Stoichiometry');
    catch ME
        getReport(ME);
        error('Failed to pass a test. Please resolve the problem before pushing.');
    end
    fprintf( 'All tests passed!\n' );
end

function doTest( dir )
    cd(dir);
    TestFeature;
    cd('..');
end