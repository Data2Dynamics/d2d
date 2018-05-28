% Load model
arInit;
arLoadModel('reduction');
arLoadData( 'out', 1, 'csv' );
arLoadData( 'out2', 1, 'csv' );
arCompileAll(true);

% Simulate and store states, sensitivities and state names
arSimu(true, false, true);
x1  = ar.model.condition(1).xExpSimu;
sx1 = ar.model.condition(1).sxExpSimu;
x2  = ar.model.condition(1).xExpSimu;
sx2 = ar.model.condition(1).sxExpSimu;
states = ar.model.x;

% Load model
arInit;
arLoadModel('reduction');

% Reduce the state space by the number of conserved moieties
arReduce;

arLoadData( 'out', 1, 'csv' );
arLoadData( 'out2', 1, 'csv' );
arCompileAll(true);

% Simulate and store states, sensitivities and state names
arSimu(true, false, true);
x1_r  = ar.model.condition(1).xExpSimu;
sx1_r = ar.model.condition(1).sxExpSimu;
sz1_r = ar.model.condition(1).szExpSimu;
x2_r  = ar.model.condition(1).xExpSimu;
sx2_r = ar.model.condition(1).sxExpSimu;
sz2_r = ar.model.condition(1).szExpSimu;
states_r = ar.model.x;
dep_r = ar.model.z;

% Assemble the sensitivities of the reduced model in such a way that the
% matrix is structured the same as in the full model.
[~,stateLoc]                        = ismember( states_r, states );
[~,depLoc]                          = ismember( dep_r, states );
sx_reconstituted1( :, stateLoc, :)   = sx1_r;
sx_reconstituted1( :, depLoc, :)     = sz1_r;
sx_reconstituted2( :, stateLoc, :)   = sx2_r;
sx_reconstituted2( :, depLoc, :)     = sz2_r;

% Show that the reduced model has the same trajectories
arPlot;

% Plot a comparison of all the sensitivities
plotComparison = 0;
if ( plotComparison )
    N = size(sx_reconstituted1, 3);
    Nx = ceil( sqrt( N ) );
    Ny = ceil( N / Nx );
    for a = 1 : N
        arSubplot(Nx, Ny, a, 'Plot1'); cla;
        plot(sx_reconstituted1(:,:,a)); hold on;
        plot(sx1(:,:,a), '--');
        plot(sx_reconstituted1(:,:,a) - sx1(:,:,a), 'k');
    end
    
    N = size(sx_reconstituted2, 3);
    Nx = ceil( sqrt( N ) );
    Ny = ceil( N / Nx );
    for a = 1 : N
        arSubplot(Nx, Ny, a, 'Plot2'); cla;
        plot(sx_reconstituted2(:,:,a)); hold on;
        plot(sx2(:,:,a), '--');
        plot(sx_reconstituted2(:,:,a) - sx2(:,:,a), 'k');
    end    
end
