% (Advanced topic) Muliple equilibration
% This file demonstrates how to equilibrate multiple conditions and use
% them as initial values for subsequent simulations
% Using this system requires some knowledge on how to find your way around
% the ar structure. See "help arSteadyState" for more information.

arInit;
arLoadModel('equilibration');
arLoadData('cond1', 1, 'csv');
arLoadData('cond2a', 1, 'csv');
arLoadData('cond2b', 1, 'csv');

% Merge plots
arMergePlotMulti(1, {'cond1', 'cond2a', 'cond2b'}, ...
                    {'condition 1', 'condition 2', 'condition 2'}, ...
                    'plot' );

% Use the event system (prerequisite for steady state sims)
ar.config.useEvents = 1;

%% Compile the model
arCompileAll(true);

arSetPars('k_basal', 0);
arSetPars('k_deg', -1);

%% Do not equilibrate prior to simulation
arSimu(true,true,true); arChi2(true);
arPlotY;

title( 'Press any key to use a single equilibration step (based on condition 1)' );
pause;

%% Equilibrate condition 1 and use that as initial value for all conditions
%  this is incorrect since condition 2 and 3 use an inhibitor
%  Note the fifth input (-1e7) which means equilibration runs from
%  t=-1*10^-7. This is to avoid simulating with the input during the steady
%  state simulation. A cleaner way to approach this would be to add a specific 
%  condition where the input is set to zero for steady state purposes.
model = 1;
sourceCondition = 1;
targetCondition = [1,2,3];
arFindInputs;
arSteadyState(ar, model, sourceCondition, targetCondition, -1e7);
arSimu(true,true,true); arChi2(true);
arPlotY;

title( 'Press any key to use proper equilibration setup (where condition 1 and 2 have their own equilibration step)' );
pause;

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents(ar); % Clears events (required!)
arFindInputs;
arSteadyState(ar, 1, 1, 1, -1e7);
arSteadyState(ar, 1, 2, [2,3], -1e7);
arSimu(true,true,true); arChi2(true);
arPlotY;