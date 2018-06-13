Setup

%%
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