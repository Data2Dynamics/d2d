function testSensis()
    global ar;
    fprintf( 2, 'INTEGRATION TEST FOR STEADY STATE BOUNDS\n' );
    
    compile = 1;
    m = 1;
    cond = 1;
    if ( compile )
        fprintf( 2, 'Loading model for steady state bounds test... ' );
        
        % Load models & data
        arInit;
        arClearCFiles;
        ar.config.checkForNegFluxes = false;
        arLoadModel('SteadyStateModel');
        arLoadData('steady',1,'csv',true);
        arCompileAll;



        ar.p(arFindPar('kd'))=-2;
        ar.p(arFindPar('kp'))=-3;
        fprintf( 2, 'PASSED\n' );
        
    end
    ar.config.rtol = 1e-10;
    ar.config.atol = 1e-10;

    fprintf( 2, 'Adding bounds and setting equilibration... ' );
    arClearEvents;
    arRemoveCustomResidual('all');
    c = 1e-8;

    %      A         B       C       D        E       F     G     H        I       J       K       L          M        N       O
    xl   = [1,       c,      c,      1,       c,      c,    1,    c,       c,      NaN,    NaN,    NaN ,      NaN,    NaN,    NaN ];
    xu   = [2,       0.05,   1,      2,       0.05,   1,    2,    0.05,    1,      NaN,    NaN,    NaN,       NaN,    NaN,    NaN ];
    logx = [c,       c,      c,      1,       1,      1,    c,    c,       c,      c,      c,      c,         c,      c,      c];

    %       J        K       L       M        N       O
    zl   = [1,       c,      c,      1,       c,      c];
    zu   = [2,       0.05,   1,      2,       0.05,   1];
    logz = [c,       c,      c,      1,       1,      1];

    arSteadyState(1, 1, 1);
    arSteadyStateBounds(m, cond, xl, xu, zl, zu, logx, logz);
    fprintf( 2, 'PASSED\n' );
    
    fprintf( 2, 'Testing correctness of sensitivities... ' );
    sresFD = qJac();
    error1 = max(max(ar.sres-sresFD)) / max(max(ar.sres));
    if (error1 < 1e-3)
        fprintf(2, '1 ');
    else
        error( 'FINAL ERROR TOO LARGE' );
    end
    
    ar.qLog10(1:2:end) = 0; ar.p(1:2:end) = 10.^ar.p(1:2:end);
    sresFD = qJac();
    error2 = max(max(ar.sres-sresFD)) / max(max(ar.sres));
    if (error2 < 1e-3)
        fprintf(2, '2 ');
    else
        error( 'FINAL ERROR TOO LARGE' );
    end
    
    ar.qLog10(1:2:end) = 1; ar.p(1:2:end) = log10(ar.p(1:2:end));
    ar.qLog10(2:2:end) = 0; ar.p(2:2:end) = 10.^ar.p(2:2:end);
    sresFD = qJac();
    error3 = max(max(ar.sres-sresFD)) / max(max(ar.sres));
    if (error3 < 1e-3)
        fprintf(2, '3 ');
    else
        error( 'FINAL ERROR TOO LARGE' );
    end    
    
    ar.qLog10(2:2:end) = 1; ar.p(2:2:end) = log10(ar.p(2:2:end));
    sresFD = qJac();
    error4 = max(max(ar.sres-sresFD)) / max(max(ar.sres));
    if (error4 < 1e-3)
        fprintf(2, '4 ... PASSED\n');
    else
        error( 'FINAL ERROR TOO LARGE' );
    end   
end

function sresFD = qJac()
    global ar;
    
    dp = 1e-8;
    arSimu(true, false, true);
    arCalcRes;
    arCollectRes(true)
    arCalcMerit(true);
    ref = ar.res + 0;
    for a = 1 : numel( ar.p )
        arPush;
        ar.p(a) = ar.p(a) + dp;
        arSimu(true, false, true);
        arCalcRes;
        arCollectRes(true)
        arCalcMerit(true);
        sresFD(:,a) = (ar.res - ref) / dp;
        arPop;
    end
    arSimu(true, false, true);
    arCalcRes;
    arCollectRes(true);
    arCalcMerit(true);
end