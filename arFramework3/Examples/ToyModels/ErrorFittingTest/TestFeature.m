fprintf( 'INTEGRATION TEST FOR ERROR PARAMETER ESTIMATION\n' );

fprintf( 2, 'LINEAR ERRORS' );
arInit;
arLoadModel('dummy');
arLoadData('nonlog', 1, 'csv');
arCompileAll;
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

fprintf( 2, 'LOGARITHMIC ERRORS' );
arInit;
arLoadModel('dummy');
arLoadData('log', 1, 'csv');
arCompileAll;
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

fprintf( 2, 'TWO-COMPONENT ERROR MODEL' );
arInit;
arLoadModel('dummy');
arLoadData('nonlog_linear', 1, 'csv');
arCompileAll;
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
else
    fprintf( 'PASSED\n' );
end

if ( (sd_rel > 0.18) || (sd_rel < 0.04) )
    error( 'ERROR NOT WITHIN TOLERANCE' );
else
    fprintf( 'PASSED\n' );
end
