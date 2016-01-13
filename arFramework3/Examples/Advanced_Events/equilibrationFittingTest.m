% (Advanced topic) Muliple equilibration
% This file demonstrates how to equilibrate multiple conditions and use
% them as initial values for subsequent simulations.
% Using this system requires some knowledge on how to find your way around
% the ar structure.

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

% Don't fit the standard deviation
ar.qFit(end)=0;

% Set the parameters to wrong values
arSetPars('k_basal', 0);
arSetPars('k_deg', -2);

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, [2,3], -1e7);

arSimu(true,true,true); arChi2(true);
arPlotY;

title('Pre-fitting (press any key to fit)');
pause;

arFit;
arSimu(true,true,true); arChi2(true);
arPlotY;
