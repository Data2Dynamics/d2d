function TestFeature()

fprintf( 'INTEGRATION TEST FOR ERROR PARAMETER ESTIMATION\n' );

fprintf( 2, 'LINEAR ERRORS... ' );
arInit;
arLoadModel('dummy');
arLoadData('nonlog', 1, 'csv');
arCompileAll(true);
ar.config.optimizer = 11;
arSetPars('init_mean',  5, 0, 0, 0, 1e3);
arSetPars('sd_est',     5, 1, 0, 0, 1e3);
arFit;
sd_abs = arGetPars('sd_est',0);

% These bounds are quite ad-hoc (only to detect big failures, not statistical tweaks)
if ( (sd_abs > 5.5) || (sd_abs < 4.5) )
    error( 'FINAL ERROR TOO LARGE' );
else
    fprintf( 'PASSED\n' );
end

fprintf( 2, 'LOGARITHMIC ERRORS... ' );
arInit;
arLoadModel('dummy');
arLoadData('log', 1, 'csv');
arCompileAll(true);
ar.config.optimizer = 11;
arSetPars('init_mean',  .25, 1, 0, 0, 1e3);
arSetPars('sd_est',     .25, 1, 0, 0, 1e3);
arFit;
sd_abs = arGetPars('sd_est',0);

% These bounds are quite ad-hoc (only to detect big failures, not statistical tweaks)
if ( (sd_abs > .26) || (sd_abs < .24) )
    error( 'FINAL ERROR TOO LARGE' );
else
    fprintf( 'PASSED\n' );
end

% Check whether qError detects the error model parameters correctly
fprintf( 2, 'CHECK WHETHER ar.qError DETECTS ERROR PARS CORRECTLY FOR ONE COMPONENT CASE... ' );
if ~(numel( intersect( ar.pLabel(ar.qError==1), {'sd_est'} ) ) == 1)
    error( 'ar.qError INCORRECTLY SPECIFIED! CHECK ARLOADMODEL/ARLOADDATA/ARLINK FOR MISTAKES' );
else
    fprintf( 'PASSED\n' );
end

fprintf( 2, 'TWO-COMPONENT ERROR MODEL... ' );
arInit;
arLoadModel('dummy');
arLoadData('nonlog_linear', 1, 'csv');
arCompileAll(true);
ar.config.optimizer = 11;
arSetPars('init_mean',   5, 0, 0, 0, 1e3);
arSetPars('sd_est_abs', 50, 1, 0, 0, 1e3);
arSetPars('sd_est_rel', .1, 1, 0, 0, 1e3);
arFit;
sd_abs = arGetPars('sd_est_abs',0);
sd_rel = arGetPars('sd_est_rel',0);

% These bounds are quite ad-hoc (only to detect big failures, not statistical tweaks)
if ( (sd_abs > 70) || (sd_abs < 30) )
    error( 'ERROR NOT WITHIN TOLERANCE' );
end

if ( (sd_rel > 0.18) || (sd_rel < 0.04) )
    error( 'ERROR NOT WITHIN TOLERANCE' );
end
fprintf( 'PASSED\n' );

% Check whether qError detects the error model parameters correctly
fprintf( 2, 'CHECK WHETHER ar.qError DETECTS ERROR PARS CORRECTLY FOR TWO COMPONENT CASE... ' );
if ~(numel( intersect( ar.pLabel(ar.qError==1), {'sd_est_abs', 'sd_est_rel'} ) ) == 2)
    error( 'FAILED!\nar.qError INCORRECTLY SPECIFIED! CHECK ARLOADMODEL/ARLOADDATA/ARLINK FOR MISTAKES' );
else
    fprintf( 'PASSED\n' );
end

fprintf( 2, 'CHECK WHETHER ar.qError DETECTS ERROR PARS CORRECTLY FOR FILENAME BASED CASE... ' );
arInit;
arLoadModel('dummyFn');
arLoadData('fn_nonlog', 1, 'csv');
arLoadData('fn_nonlog_linear', 1, 'csv');
arCompileAll(true);

% Check whether qError detects the error model parameters correctly
if ( ~(numel( intersect( ar.pLabel(ar.qError==1), {'sd_est_abs_fn_nonlog', 'sd_est_abs_fn_nonlog_linear', 'sd_est_rel_fn_nonlog', 'sd_est_rel_fn_nonlog_linear'} ) ) == 4) ) || (sum(ar.qError==1) ~= 4)
    error( 'FAILED!\nar.qError INCORRECTLY SPECIFIED! CHECK ARLOADMODEL/ARLOADDATA/ARLINK FOR MISTAKES' );
else
    fprintf( 'PASSED\n' );
end

fprintf( 2, 'CHECK WHETHER D2D CORRECTLY EXCLUDES DYNAMIC AND OBSERVATION PARAMETERS FROM ar.qError... ' );
arInit;
arLoadModel('withscale');
arLoadData('fn_nonlog', 1, 'csv');
arLoadData('fn_nonlog_linear', 1, 'csv');
arCompileAll(true);

% Check whether qError detects the error model parameters correctly
if ~(numel( intersect( ar.pLabel(ar.qError==1), {'sd_est_abs_fn_nonlog', 'sd_est_abs_fn_nonlog_linear', 'sd_est_rel_fn_nonlog', 'sd_est_rel_fn_nonlog_linear'} ) ) == 4 ) || (sum(ar.qError==1) ~= 4)
    error( 'FAILED!\nar.qError INCORRECTLY SPECIFIED! CHECK ARLOADMODEL/ARLOADDATA/ARLINK FOR MISTAKES' );
else
    fprintf( 'PASSED\n' );
end

fprintf( 2, 'CHECK WHETHER D2D CORRECTLY HANDLES RANDOMS WITH RESPECT TO ar.qError ... ' );
arInit;
arLoadModel('nExpID');
arLoadData('fn_log', 1, 'csv');
arCompileAll(true);

% Check whether qError detects the error model parameters correctly
if ~(numel( intersect( ar.pLabel(ar.qError==1), {'sd_est_abs_fn_log_nExpID1', 'sd_est_abs_fn_log_nExpID2', 'sd_est_rel_fn_log_nExpID1', 'sd_est_rel_fn_log_nExpID2'} ) ) == 4 ) || (sum(ar.qError==1) ~= 4)
    error( 'FAILED!\nar.qError INCORRECTLY SPECIFIED! CHECK ARLOADMODEL/ARLOADDATA/ARLINK FOR MISTAKES' );
else
    fprintf( 'PASSED\n' );
end

