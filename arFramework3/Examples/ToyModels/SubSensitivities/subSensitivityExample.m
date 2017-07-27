% Subsensitivity

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

% Tight tolerances are required for this system
ar.config.rtol=1e-9;
ar.config.atol=1e-10;

% Simulate and show
arSimu(true,true,true); arGetMerit;
arPlotY;

% Only fit a subset of the parameters
ar.qFit=ones(size(ar.qFit)); ar.qFit([2,4,6,8,10,12,14,16])=0;

ar.config.sensitivitySubset=0;
tic;
for a = 1 : 15
    arSimu(true, true, true);
end
arCalcMerit(true);
sres_without = ar.sres;
q = toc;
fprintf( 'Time without subsensitivities: %g\n', q );

ar.config.sensitivitySubset=1;
tic;
for a = 1 : 15
    arSimu(true, true, true);
end
arCalcMerit(true);
q = toc;
sres_with = ar.sres;
fprintf( 'Time with subsensitivities: %g\n', q );

% Verify that sres is the same
fprintf( 'Sum of squares difference in sensitivities %g\n', sum( sum( (sres_with(:, ar.qFit==1) - sres_without(:, ar.qFit==1)).^2 ) ) );
