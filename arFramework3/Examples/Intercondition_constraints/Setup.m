% This example file demonstrates how to regularize predictions in such a
% way as to minimize differences between different conditions

arInit;
arLoadModel( 'min_model' );
arLoadData( 'drug', 1, 'csv' );
arLoadData( 'noDrug', 1, 'csv'  );
arLoadData( arGenerateConditionData( 'steadystate', 'input_drug', 0, 'input_stimulus', 0 ), 1, 'csv' );
arCompileAll;

C = arFindCondition('', 'input_drug', '0', 'input_stimulus', '0' );
arSteadyState(1, C, 'all');

% Initial parameters
ar.p = -1 * ones(size(ar.p));
ar.p(arFindPar('kq')) = 3;      % Deliberately set this one wrong

%% Fit and plot without constraints
%  Note that the two states are estimated as being different, while they don't have to be
arFit;
arSimu(false, true); arChi2;
ar.model(1).qPlotXs = [1 1 0];
figure; subplot(2,1,1);
C1 = arFindCondition('', 'input_drug', '0', 'input_stimulus', '1');
C2 = arFindCondition('', 'input_drug', '10', 'input_stimulus', '1');
s  = find(ismember( ar.model.x, 'C_state' ));
plot( ar.model(1).condition(C1).tFine, ar.model(1).condition(C1).xFineSimu( :, s ), 'k', 'LineWidth', 2 );
hold on;
plot( ar.model(1).condition(C2).tFine, ar.model(1).condition(C2).xFineSimu( :, s ), 'r', 'LineWidth', 2 );
hold on;
plot( ar.model(1).condition(C2).tFine, ar.model(1).condition(C1).xFineSimu( :, s )-ar.model(1).condition(C2).xFineSimu( :, s ), 'k--', 'LineWidth', 2 );
title( sprintf( 'The data does not require these two to be different.\nPress a key to introduce condition constraint.' ) );
legend( 'State3', 'State4', 'Difference' );

pause;

%% Add soft relative condition constraints
% Fetch relevant conditions
m1 = 1;
m2 = 1;
C1 = arFindCondition('', 'input_drug', '0', 'input_stimulus', '1');
C2 = arFindCondition('', 'input_drug', '10', 'input_stimulus', '1');

% Select time points to use for constraints
t = 0:1:120;

% Sigma per time point (relative units; extremely soft)
sigma = 1;

% States to constrain (assume we don't know which ones to constrain)
states = 'all';

%% Add it
arAddConditionConstraint( m1, C1, m2, C2, t, sigma, states );

%% Fit and plot
arFit; arSimu(false, true); arChi2;
subplot(2,1,2);
C1 = arFindCondition('', 'input_drug', '0', 'input_stimulus', '1');
C2 = arFindCondition('', 'input_drug', '10', 'input_stimulus', '1');
s  = find(ismember( ar.model.x, 'C_state' ));
plot( ar.model(1).condition(C1).tFine, ar.model(1).condition(C1).xFineSimu( :, s ), 'k', 'LineWidth', 2 );
hold on;
plot( ar.model(1).condition(C2).tFine, ar.model(1).condition(C2).xFineSimu( :, s ), 'r', 'LineWidth', 2 );
hold on;
plot( ar.model(1).condition(C2).tFine, ar.model(1).condition(C1).xFineSimu( :, s )-ar.model(1).condition(C2).xFineSimu( :, s ), 'k--', 'LineWidth', 2 );
title( 'With condition constraint between the states' );
legend( 'State3', 'State4', 'Difference' );

