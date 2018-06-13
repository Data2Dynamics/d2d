% Example demonstrating the use of the fast sensitivity system

arInit;

% This flag makes sure that the model is compiled with root finding
% capabilities.
ar.config.fastEquilibration = 1;

arLoadModel('equilibration2');

% If you want to use fast equilibration, this step is important if you have 
% conserved moieties! This is to ensure that dfdx is full rank.
arReduce;

% Load the data
arLoadData( 'reduced_ss_condi1', 1, 'csv' );
arLoadData( 'reduced_ss_condi2', 1, 'csv' );
arLoadData( 'reduced_ss_condi3', 1, 'csv' );

% Merge plots
arMergePlotMulti(1, {'reduced_ss_condi1', 'reduced_ss_condi2', 'reduced_ss_condi3'}, ...
                    {'condition 1', 'condition 2', 'condition 2'}, ...
                     'plot' );

% Use the event system (prerequisite for steady state sims)
ar.config.useEvents = 1;

%% Compile the model
arCompileAll(true);

