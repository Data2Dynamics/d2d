function TestFeature()

fprintf( 2, 'INTEGRATION TEST FOR SUBSENSITIVITIES\n' );

fprintf( 2, 'Loading model for subsensitivity test... ' );
% Subsensitivity test
arInit;
arLoadModel('subsensi');
arLoadData('wt', 1, 'csv');
arLoadData('condi1', 1, 'csv');
arLoadData('condi2', 1, 'csv');

% Merge plots
arMergePlotMulti(1, {'wt', 'condi1', 'condi2'}, ...
                    {'wt', 'Lower degradation', 'Lower basal rate'}, ...
                    'plot' );

% Use the event system (prerequisite for steady state sims)
ar.config.useEvents = 1;
ar.config.maxsteps = 1e5;

%% Compile the model
arCompileAll(true);

% True parameters
ar.p( arFindPar( 'amount' ) ) = -1;
ar.p( arFindPar( 'hill' ) ) = 1;
ar.p( arFindPar( 'init_stateB' ) ) = -1;
ar.p( arFindPar( 'input_bas_inh' ) ) = -1;
ar.p( arFindPar( 'k_C' ) ) = 2;
ar.p( arFindPar( 'k_D' ) ) = -1;
ar.p( arFindPar( 'k_b' ) ) = -1;
ar.p( arFindPar( 'k_basal' ) ) = -1;
ar.p( arFindPar( 'k_degA' ) ) = -1;
ar.p( arFindPar( 'k_degB' ) ) = -0.7;
ar.p( arFindPar( 'k_degC' ) ) = -0.5;
ar.p( arFindPar( 'k_degD' ) ) = -0.5;
ar.p( arFindPar( 'k_inh' ) ) = 3;
ar.p( arFindPar( 'ka_C' ) ) = 1;
ar.p( arFindPar( 'sd_' ) ) = -3;

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents; % Clears events (required!)
arFindInputs;
arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, 2, -1e7);
arSimu(true,true,true); arGetMerit;
ar.config.rtol=1e-9;
ar.config.atol=1e-10;
fprintf( 'PASSED\n' );


fprintf( 2, 'Testing subset patterns... ' );

patterns{1} = [2,4,6,8,10,12,14,16];
patterns{2} = 1:16;
patterns{3} = [2,3,6,7,8,9,10,12,14,16];
patterns{4} = [12,13,14];
patterns{5} = [1,3,4,7];
patterns{6} = [5,9,11];
patterns{7} = [1,2,3];
patterns{8} = [2,3,4];
patterns{9} = [6,7,9];

for a = 1 : numel( patterns )
    % Only fit a subset of the parameters
    ar.qFit=ones(size(ar.qFit));
    ar.qFit(patterns{a})=0;

    ar.config.sensitivitySubset=0;
    arSimu(true, true, true); arCalcMerit(true);
    sres_without = ar.sres + 0;

    ar.config.sensitivitySubset=1;
    arSimu(true, true, true); arCalcMerit(true);
    sres_with = ar.sres + 0;
    
    fprintf( 2, '%d ', a );
    
    % Verify that sres is the same
    diff = sres_with(:,ar.qFit==1) - sres_without(:,ar.qFit==1);
    if ( sum( sum( (diff).^2 ) ) > ar.config.atol * 1000 )
        error( 'FAILED AT PATTERN %d! Error was: %d', a, sum( sum( (diff).^2 ) ) );
    end
end

fprintf( 2, 'PASSED\n' );