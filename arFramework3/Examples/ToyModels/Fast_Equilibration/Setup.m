% Example demonstrating the use of the fast sensitivity system

arInit;

arLoadModel('equilibration2');

% This step is important if you have conserved moieties!
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

% Enable fast sensitivities
arFastSensis;

% Don't fit the standard deviation
ar.qFit(end)=0;
ar.p(6)=0

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, [2,3], -1e7);

% Compare performance
comparePerformance = 0;
ar.config.rtol=1e-9; 
ar.config.atol=1e-9;
if comparePerformance
    ar.config.turboSSSensi=0;
    tic
        % Set the parameters to wrong values
        ar.p = -ones(size(ar.p));
        arCalcMerit;
        arFit;
    toc

    ar.config.turboSSSensi=1;
    tic
        % Set the parameters to wrong values
        ar.p = -ones(size(ar.p));
        arCalcMerit;
        arFit;
    toc
end