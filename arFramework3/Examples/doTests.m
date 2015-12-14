% Please verify any changes made to arCompile, arCompileAll, arLoadModel,
% arLoadData, arSimu, arSimuCalc, arFit by running these integration tests.
%
% Usage: doTests( list of tests )
%
% Examples:
%       doTests                                             - Performs all tests
%       doTests('Volume_Fitting_Test', 'Splines')           - Perform specified two tests
%
function doTests( varargin )
    global ar;
    global arOutputLevel;

    fprintf(2, 'This collection of tests checks whether specific functions\n');
    fprintf(2, 'in the D2D functions are working correctly. Please run this\n');
    fprintf(2, 'routine before pushing any changes to internal functions in D2D\n' );
    fprintf(2, 'to reduce the risk of pushing code that breaks existing\nfunctionality.\n\n' );

    tests = {   'Advanced_Events', 'Volume_Estimation', 'Splines', ...
                'Stoichiometry', 'Step_Estimation' };
    
    if ( nargin > 0 )
        activeTests = argSwitch( tests, varargin{:} );
    else
        activeTests = argSwitch( tests, tests{:} );
    end
            
    try
        % These tests require lsqnonlin
        if ( license('checkout', 'Optimization_Toolbox') )
            ar.config.optimizer = 1;
            fprintf( 'Optimization toolbox found, using lsqnonlin.\n');
        else
            ar.config.optimizer = 4;
            fprintf( 'Optimization toolbox not found, switching optimizer to STRSCNE.\n');
        end
        
        for a = 1 : length( tests )
            if ( activeTests.(tests{a}) )
                doTest( tests{a} );
            else
                fprintf(2, 'SKIPPING %s\n', tests{a} );
            end
        end
    
        fprintf('Testing complete!\n');
        arOutputLevel = 1;
    catch
        fprintf('Testing failed: please resolve errors and try again.\n');
        arOutputLevel = 1;
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

function [opts] = argSwitch( switches, varargin )

    for a = 1 : length(switches)
        opts.(switches{a}) = 0;
    end

    a = 1;
    while (a <= length(varargin))
        if ( max( strcmpi( varargin{a}, switches ) ) == 0 )
            str = sprintf( 'Legal switch arguments are:\n' );
            str = [str sprintf( '%s\n', switches{:} ) ];%#ok<AGROW>
            error( 'Invalid switch argument was provided. Provided %s, %s', varargin{a}, str );
        else
            fieldname = switches{ find( strcmpi( varargin{a}, switches ) ) };
        end
        
        val = 1;
        if ( length(varargin) > a )
            if isnumeric( varargin{a+1} )
                val = varargin{a+1};
                a = a + 1;
            end
        end
        
        opts.(fieldname) = val;
        a = a + 1;
    end
end