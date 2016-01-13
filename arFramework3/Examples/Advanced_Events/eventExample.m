% This file demonstrates construction of a model with complex events as 
% well as pre-equilibration of a condition.

arInit;
arLoadModel('events');
ar.config.useEvents = 1;

%% Compile the model
arCompileAll(true);

arSetPars('k_prod', .2);
arSetPars('k_deg', -3);
arSetPars('k_tr', -0.5);

%% Add some on the fly event state changes (no recompile necessary)
%  State changes are of the form aX+b where X corresponds to the state
%  variable value before the change.
arAddEvent(1, 1, 50,  'stateA', 1,  10);    % Add 10A at t=50 [A=(1*A+10)]
arAddEvent(1, 1, 60,  'stateA', 0.5 );      % Double the volume (half the concentration) at t=60 (A=(0.5*A))
arAddEvent(1, 1, 60,  'stateB', 0.5 );      % Doubling the volume also has to be taken into account for state B (B=(0.5*B))
arAddEvent(1, 1, 70,  'stateA', 1,  10);    % Add some A at t=70
arAddEvent(1, 1, 80,  'stateA', 0.5 );      % Double the volume (half the concentration) at t=80
arAddEvent(1, 1, 80,  'stateB', 0.5 );
arAddEvent(1, 1, 90,  'stateA', 1,  10);    % Add some A at t=90
arAddEvent(1, 1, 100, 'stateA', 0.5 );      % Double the volume (half the concentration) at t=100
arAddEvent(1, 1, 100, 'stateB', 0.5 );
arAddEvent(1, 1, 110, 'stateA', 1,  10);    % Add some A at t=110
arAddEvent(1, 1, 120, 'stateA', 0.5 );      % Double the volume (half the concentration) at t=110
arAddEvent(1, 1, 120, 'stateB', 0.5 );

arSimu(true,true,true); arChi2(true);
arPlotX;
title( 'Press any key to show pre-equilibrated simulation' );
pause;

%% Equilibrate condition 1 and use that as initial value for condition 1
model = 1;
sourceCondition = 1;
targetCondition = 1;
arSteadyState(model, sourceCondition, targetCondition);

arSimu(true,true,true); arChi2(true);
arPlotX;