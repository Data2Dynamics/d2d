% (Advanced topic) Equilibration using rootfinding
arInit;
arLoadModel('conservedPools');
arLoadData('steadystate');
arLoadData('nonsteadystate');

%% Compile the model
arCompileAll(true);

% Set step amount
ar.p(3) = 2;
arSteadyState(1,arFindCondition('steadystate', 'exact'), arFindCondition('nonsteadystate', 'exact'));

%% Do not equilibrate prior to simulation
arSimu(true,true,true); arChi2(true);

% The rootfinding should conserve the totals
xRoot = arFindRoots;
volRatio = (.2/4);
total_S1 = xRoot( ismember(ar.model.x, 'S1') ) + xRoot( ismember(ar.model.x, 'pS1') );
total_S2 = xRoot( ismember(ar.model.x, 'S2') ) + xRoot( ismember(ar.model.x, 'pS2') ) + 2 * xRoot( ismember(ar.model.x, 'pS2_pS2') ) + xRoot( ismember(ar.model.x, 'ppS2') ) * volRatio;
fprintf( 'S1 %g (should be 5)\nS2 %g (should be 8)\n\n', total_S1, total_S2 );

% Deliberately make it misfit
ar.config.rootfinding = 0;
ar.p=.01*ones(size(ar.p));
tic;
arFit;
sstime = toc;
fprintf( 'Time taken with simulation %g\n', sstime );

% Deliberately make it misfit
ar.config.rootfinding = 1;
ar.p=.01*ones(size(ar.p));
tic;
arFit;
rftime = toc;
fprintf( 'Time taken with rootfinding %g\n', rftime );