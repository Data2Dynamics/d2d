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

% Enable fast sensitivities
arFastSensis;

% Don't fit the standard deviation
ar.qFit(end)=0;
ar.p(6)=0;

%% Equilibrate condition 1 and use that as initial value for condition 1
%  Equilibrate condition 2 and use that as initial condition for 2 and 3
arClearEvents; % Clears events
arFindInputs;
arSteadyState(1, 1, 1, -1e7);
arSteadyState(1, 2, [2,3], -1e7);

% Compare performance
comparePerformance = 1;
ar.config.rtol=1e-9; 
ar.config.atol=1e-9;
if comparePerformance
    % This is the default 'slow' equilibration procedure
    ar.config.rootfinding = 0;
    ar.config.turboSSSensi = 0;
    disp( 'Normal equilbration' );
    tic
        % Set the parameters to wrong values
        ar.p = -ones(size(ar.p));
        arCalcMerit;
        arFit;
    toc

    % This simulates the system but computes the steady state sensitivities
    % implicitly (requires full rank dfdx i.e. no conserved moieties)
    ar.config.rootfinding = 0;
    ar.config.turboSSSensi = 1;
    disp( 'Simulate core ODEs only for equilibration, implicitly calculate sensitivities' );
    tic
        % Set the parameters to wrong values
        ar.p = -ones(size(ar.p));
        arCalcMerit;
        arFit;
    toc
end