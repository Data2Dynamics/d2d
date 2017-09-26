function TestFeature()

global ar;

fprintf( 2, 'INTEGRATION TEST FOR FAST EQUILIBRATION\n' );

fprintf( 2, 'Loading model for fast equilibration test (not reduced) ... ' );
arInit;
ar.config.fastEquilibration = 0;
arLoadModel('equilibration2');
fprintf( 2, ' [ OK ]\n' );

fprintf( 2, 'Testing whether it correctly throws an error when attempting to use fast equilibration ...' );

% This should throw!
thrown = 0;
try
    arFastSensis
catch
    thrown = 1;
end
if ( thrown )
    fprintf( 2, '[ OK ]\n' );
else
    error( 'Attempt to use fast steady state sensis should have thrown an error on a model with conserved moieties' );
end

fprintf( 2, 'Loading model for fast equilibration test (with data, reduced) ... ' );
arInit;
ar.config.fastEquilibration = 0;
arLoadModel('equilibration2');
arReduce;
arLoadData( 'reduced_ss_condi1', 1, 'csv' );
arLoadData( 'reduced_ss_condi2', 1, 'csv' );
arLoadData( 'reduced_ss_condi3', 1, 'csv' );
arCompileAll(true);
fprintf( 2, '[ OK ]\n' );

fprintf( 2, 'Activating fast sensitivities ... ' );
arFastSensis
fprintf( 2, '[ OK ]\n' );

arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, [2,3], -1e7);    

fprintf( 2, 'Simulating sensitivities implicitly and explicitly ... ' );
ar.config.rtol = 1e-10; ar.config.atol = 1e-10;
rtol = 1e-8; atol = 1e-8;
for d = 1 : numel( ar.model.ss_condition )
    ar.config.turboSSSensi=1;
    arCalcMerit;
    sxEq = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;

    ar.config.turboSSSensi=0;
    arCalcMerit;
    subplot(numel(ar.model.ss_condition), 1, d);
    sxTrue = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;        
end

fail = sum( (((sxEq-sxTrue)./(sxTrue))>rtol) & ((sxEq-sxTrue)>atol) );
if ( fail )
    error( 'Failure: Sensitivity difference is too large' );
else
    fprintf( 2, '[ OK ]\n' );
end

fprintf( 2, 'Loading model for fast equilibration test (no data, reduced) ... ' );
arInit;
ar.config.fastEquilibration = 0;
arLoadModel('equilibration2');
arReduce;
arCompileAll(true);

fprintf( 2, 'Activating fast sensitivities ... ' );
arFastSensis
fprintf( 2, '[ OK ]\n' );

arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);  

fprintf( 2, 'Simulating sensitivities implicitly and explicitly ... ' );
ar.config.rtol = 1e-10; ar.config.atol = 1e-10;
rtol = 1e-8; atol = 1e-8;
for d = 1 : numel( ar.model.ss_condition )
    ar.config.turboSSSensi=1;
    arCalcMerit;
    sxEq = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;

    ar.config.turboSSSensi=0;
    arCalcMerit;
    subplot(numel(ar.model.ss_condition), 1, d);
    sxTrue = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;        
end

fail = sum( (((sxEq-sxTrue)./(sxTrue))>rtol) & ((sxEq-sxTrue)>atol) );
if ( fail )
    error( 'Failure: Sensitivity difference is too large' );
else
    fprintf( 2, '[ OK ]\n' );
end

fprintf( 2, 'Loading model for fast equilibration test with rootfinding in C (no data, reduced) ... ' );
arInit;
ar.config.fastEquilibration = 1;
arLoadModel('equilibration2');
arReduce;
arCompileAll(true);

fprintf( 2, 'Activating fast sensitivities ... ' );
arFastSensis
fprintf( 2, '[ OK ]\n' );

arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);  

fprintf( 2, 'Simulating sensitivities implicitly and explicitly ...\n' );
ar.config.rtol = 1e-10; ar.config.atol = 1e-10;
rtol = 1e-8; atol = 1e-8;
for d = 1 : numel( ar.model.ss_condition )
    ar.config.C_rootfinding=0;
    ar.config.rootfinding=1;
    arCalcMerit;
    sxEq = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;
    
    ar.config.C_rootfinding=1;
    ar.config.rootfinding=1;
    arCalcMerit;
    sxEqC = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;    

    ar.config.rootfinding=0;
    arCalcMerit;
    subplot(numel(ar.model.ss_condition), 1, d);
    sxTrue = reshape(ar.model.ss_condition(d).sxFineSimu(end,:,:), 1, numel(ar.model.ss_condition(d).sxFineSimu(end,:,:))) + 0;        
end

fprintf( 2, 'Testing whether sensitivities of solution obtained with rootfinding in MATLAB are identical to those with simulation ... ' );
fail = sum( (((sxEq-sxTrue)./(sxTrue))>rtol) & ((sxEq-sxTrue)>atol) );
if ( fail )
    error( 'Failure: Sensitivity difference is too large between simulated and rootfinding in MATLAB' );
else
    fprintf( 2, '[ OK ]\n' );
end

fprintf( 2, 'Testing whether sensitivities of solution obtained with rootfinding in C are identical to those with simulation ... ' );
fail = sum( (((sxEqC-sxTrue)./(sxTrue))>rtol) & ((sxEqC-sxTrue)>atol) );
if ( fail )
    error( 'Failure: Sensitivity difference is too large between simulated and rootfinding in C' );
else
    fprintf( 2, '[ OK ]\n' );
end
