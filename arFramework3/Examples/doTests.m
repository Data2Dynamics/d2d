% Please verify any changes made to arCompile, arCompileAll, arLoadModel,
% arLoadData, arSimu, arSimuCalc, arFit by running these integration tests.
%
% Usage: doTests( list of tests )
%
% Examples:
%       doTests                                             - Performs all tests
%       doTests('Volume_Fitting_Test', 'Splines')           - Perform specified two tests
%
% If you also wish to run the long tests, call doTests('long')
%
% The available tests are:
%     'SubSensitivities', 'Advanced_Events', 'Volume_Estimation', 'Splines', 'DataFilterTest',
%     'Stoichiometry', 'DallaMan2007_GlucoseInsulinSystem', 'Step_Estimation', 
%     'ErrorFittingTest', 'Flux_Estimation', 'MultiCondition_Test', 'TurboSplines', 
%     'ResponseCurve', 'PreProcessorTest', 'State_Reduction', 'Fast_Equilibration',  
%     'Predictor_Test', 'FieldTester', 'SteadyStateBounds', 'DataFilterTest', 
%     'Benchmark_Simu_Test', 'Boehm_JProteomeRes2014', 'PEtab_Tests'
function doTests( varargin )
    global ar;
    global arOutputLevel;
    global arStrict;
    
    arOutputLevel = 0;
    arStrict      = 1;
    warning off;

    fprintf(2, 'This collection of tests checks whether specific functions\n');
    fprintf(2, 'in the D2D functions are working correctly. Please run this\n');
    fprintf(2, 'routine before pushing any changes to internal functions in D2D\n' );
    fprintf(2, 'to reduce the risk of pushing code that breaks existing\nfunctionality.\n\n' );

    tests = {   'SubSensitivities', 'Advanced_Events', 'Volume_Estimation', 'Splines', 'DataFilterTest', ...
                'Stoichiometry', 'DallaMan2007_GlucoseInsulinSystem', 'Step_Estimation', ...
                'ErrorFittingTest', 'Flux_Estimation', 'MultiCondition_Test', 'TurboSplines', ...
                'ResponseCurve', 'PreProcessorTest', 'State_Reduction', 'Fast_Equilibration', ... 
                'Predictor_Test', 'FieldTester', 'SteadyStateBounds', 'SD_Test', ...
                'Benchmark_Simu_Test', 'Boehm_JProteomeRes2014', 'PEtab_Tests'};
            
    longtests = { 'Benchmark_Simu_Test' };
    
    dependencies = { {}, {}, {}, {}, {}, {}, {'TranslateSBML'}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {'TranslateSBML'}, {'TranslateSBML'}};
    
    if ( nargin > 0 && strcmp( varargin{1}, 'long' ) )
        varargin = setdiff( varargin, 'long' );
        activeTests = argSwitch( tests, tests{:} );
    else
        shortTests = setdiff( tests, longtests );
        activeTests = argSwitch( tests, shortTests{:} );
    end    
    
    if ( numel(varargin) > 0 )
        activeTests = argSwitch( tests, varargin{:} );
    end
    activeTests
    successes = 0;
    failures = 0;
    skipped = 0;
    fails = zeros( numel( tests ), 1 );
    
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
            dependenciesFound = 1;
            for jd = 1 : length( dependencies{a} )
                if ( exist(dependencies{a}{jd}, 'file') == 0 )
                    dependenciesFound = 0;
                    fprintf('SKIPPING %s (missing dependency %s)\n', tests{a}, dependencies{a}{jd});
                    skipped = skipped + 1;
                end
            end
            if ( dependenciesFound )
                try
                    doTest( tests{a} );
                    successes = successes + 1;
                catch
                    fails(a) = 1;
                    failures = failures + 1;
                    fprintf( '\n\n' );
                end 
            end
        else
            fprintf(2, 'SKIPPING %s\n', tests{a} );
        end
    end
   
    fprintf('\n\n----------------\nTesting complete! %d test%s passed, %d test%s failed, %d test%s skipped due to missing dependencies.\n', successes, pluralize(successes), failures, pluralize(failures), skipped, pluralize(skipped) );
    if ( failures > 0 )
        str = sprintf( '%s ', tests{fails==1} );
        fprintf( 2, 'The following tests failed: %s\n\n********************************************************************************\n* PLEASE DO *NOT* COMMIT your changes to the master branch before fixing your  *\n* changes in such a way that all tests succeed.                                *\n********************************************************************************\n\n', str );
    end
    
    arOutputLevel = 2;
    arStrict = 0;
    warning on;
end

function r = pluralize(value)
    if (value~=1)
        r = 's';
    else
        r = '';
    end
end

function doTest( dir )
    d = pwd;
    example_folder = strrep(which('arInit'),'arInit.m','Examples');
    if isunix
        command = sprintf('find %s -type d -name "%s"',example_folder,dir);
    else
        command = sprintf('dir "%s" /AD /b /s | findstr /r "\\\\%s$"',example_folder,dir);
    end
    [~,dir] = system(command);
    dir = strread(dir, '%s', 'delimiter', sprintf('\n'));
    cd(dir{1});
    
    % Suppress output
    global arOutputLevel;
    arOutputLevel = 0;
    try
    TestFeature; cd(d);
    arClearCFiles;
    drawnow;
    sprintf('\n');
    
    catch
        try
            fprintf(2, 'UNIT TEST FAILED: RESTARTING WITH VERBOSE OUTPUT\n\n' );
            arOutputLevel = 2;
            TestFeature;
        catch ME
            fprintf(getReport(ME));
            arOutputLevel = 0;
            cd(d);
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
            fieldname = switches{ strcmpi( varargin{a}, switches ) };
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