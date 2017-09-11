function TestFeature()

fprintf( 2, 'INTEGRATION TEST FOR AUTOMATIC MODEL STATE REDUCTION\n' );

fprintf( 2, 'Loading and compiling model (no reduction) ... ' );
% Load model
arInit;
arLoadModel('reduction');
arLoadData( 'out', 1, 'csv' );
arLoadData( 'out2', 1, 'csv' );
arCompileAll(true);
fprintf( 2, 'PASSED\n' );

fprintf( 2, 'Simulating model (no reduction) ... ' );
% Simulate and store states, sensitivities and state names
arSimu(true, false, true);
x1  = ar.model.condition(1).xExpSimu;
sx1 = ar.model.condition(1).sxExpSimu;
x2  = ar.model.condition(2).xExpSimu;
sx2 = ar.model.condition(2).sxExpSimu;
states = ar.model.x;
fprintf( 2, 'PASSED\n' );

fprintf( 2, 'Loading model (reduction) ... ' );
% Load model
arInit;
arLoadModel('reduction');
fprintf( 2, 'PASSED\n' );

% Reduce the state space by the number of conserved moieties
fprintf( 2, 'Reducing model ... ' );
arReduce;
fprintf( 2, 'PASSED\n' );

fprintf( 2, 'Compiling reduced model ... ' );
arLoadData( 'out', 1, 'csv' );
arLoadData( 'out2', 1, 'csv' );
arCompileAll(true);
fprintf( 2, 'PASSED\n' );

% Simulate and store states, sensitivities and state names
fprintf( 2, 'Simulating model (reduction) ... ' );
arSimu(true, false, true);
x1_r  = ar.model.condition(1).xExpSimu;
sx1_r = ar.model.condition(1).sxExpSimu;
sz1_r = ar.model.condition(1).szExpSimu;
x2_r  = ar.model.condition(2).xExpSimu;
sx2_r = ar.model.condition(2).sxExpSimu;
sz2_r = ar.model.condition(2).szExpSimu;
states_r = ar.model.x;
dep_r = ar.model.z;
fprintf( 2, 'PASSED\n' );

% Tolerances
reltol = 1e-6;
abstol = 1e-4;

fprintf( 2, 'Evaluating data description of reduced model, main model condition...' );
if ( (sum(sum((ar.model.data(1).yExp-ar.model.data(1).yExpSimu).^2)) < abstol ) )
    fprintf( 2, 'PASSED\n' );
else
    error( 'FAILED: Model error too large!\n' );
end

fprintf( 2, 'Evaluating data description of reduced model, derived model condition...' );
if ( (sum(sum((ar.model.data(2).yExp-ar.model.data(2).yExpSimu).^2)) < abstol ) )
    fprintf( 2, 'PASSED\n' );
else
    error( 'FAILED: Model error too large!\n' );
end

% Assemble the sensitivities of the reduced model in such a way that the
% matrix is structured the same as in the full model.
fprintf( 2, 'Comparing sensitivities of reduced model to full model ...' );
[~,stateLoc]                        = ismember( states_r, states );
[~,depLoc]                          = ismember( dep_r, states );
sx_reconstituted1                   = [];
sx_reconstituted2                   = [];
sx_reconstituted1( :, stateLoc, :)  = sx1_r;
sx_reconstituted1( :, depLoc, :)    = sz1_r;
sx_reconstituted2( :, stateLoc, :)  = sx2_r;
sx_reconstituted2( :, depLoc, :)    = sz2_r;

% Verify that the sensitivities are the same
for ji = 1 : size( sx1, 1 )
    for jj = 1 : size( sx1, 2 )
        for jk = 1 : size( sx1, 3 )
            exceed1(ji, jj, jk) = min( [ sx_reconstituted1(ji, jj, jk) - sx1(ji, jj, jk) > abstol, ( ( sx_reconstituted1(ji, jj, jk) - sx1(ji, jj, jk) ) / sx1(ji, jj, jk) ) > reltol ] );
        end
    end
end

% Verify that the sensitivities are the same
for ji = 1 : size( sx2, 1 )
    for jj = 1 : size( sx2, 2 )
        for jk = 1 : size( sx2, 3 )
            exceed2(ji, jj, jk) = min( [ sx_reconstituted2(ji, jj, jk) - sx2(ji, jj, jk) > abstol, ( ( sx_reconstituted2(ji, jj, jk) - sx2(ji, jj, jk) ) / sx2(ji, jj, jk) ) > reltol ] );
        end
    end
end

if ( max( max( max( exceed1 ) ) ) == 0 )
    fprintf( 2, '[ Main ]' );
    if ( max( max( max( exceed2 ) ) ) == 0 )
        fprintf( 2, ' [ Derived ] ... PASSED\n' );
    else
        error( 'FAILED: Sensitivity error too large in derived model condition!\n' );
    end
else
    error( 'FAILED: Sensitivity error too large in main model condition!\n' );
end

