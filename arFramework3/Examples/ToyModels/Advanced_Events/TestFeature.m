function TestFeature()

global ar;

fprintf( 2, 'INTEGRATION TEST FOR EQUILIBRATION\n' );

fprintf( 2, 'Loading model for equilibration test... ' );
arInit;
arLoadModel('equilibration');
arLoadData('cond1', 1, 'csv');
arLoadData('cond2a', 1, 'csv');
arLoadData('cond2b', 1, 'csv');

% Use the event system (prerequisite for steady state sims)
ar.config.useEvents = 1;

%% Compile the model
arCompileAll(true);

% Don't fit the standard deviation
ar.qFit(end)=0;

% Set the parameters to wrong values
arSetPars('k_basal', 0);
arSetPars('k_deg', -2);
fprintf( 'PASSED\n' );

fprintf( 2, 'Setting equilibration events... ' );

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents; % Clears events
arFindInputs;

arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, [2,3], -1e7);
fprintf( 'PASSED\n' );

arFit;
fprintf( 2, 'Testing fitting with equilibration event... ' );
if ((norm(ar.model.data(1).res)+norm(ar.model.data(2).res)+norm(ar.model.data(3).res))<0.01)
    fprintf(2, 'PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end

arSetPars('k_basal', 0);
arSetPars('k_deg', -1);
arSimu(false,false,true);
fprintf( 2, 'Testing correct parameters with equilibration event... ' );
if ((norm(ar.model.data(1).res)+norm(ar.model.data(2).res)+norm(ar.model.data(3).res))<0.01)
    fprintf(2, 'PASSED\n');
else
    error( 'FINAL ERROR TOO LARGE' );
end