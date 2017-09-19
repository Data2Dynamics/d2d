%
%   arFastSensis
%
%   Enables or disables fast sensitivity computation. For the computation
%   of fast sensitivities, the RHS must be full rank. For this, the model
%   has to be reduced to a form where there are no more conserved moieties.
%   Use arReduce after model loading, but before data loading to accomplish
%   this.
%

function arFastSensis()

    global ar;
    
    if ( ~isfield( ar.config, 'turboSSSensi' ) )
        ar.config.turboSSSensi = 0;
    end
    
    if ( ar.config.turboSSSensi == 1 )
        arFprintf(2, 'Disabled fast sensitivity computation\n')
        ar.config.turboSSSensi = 0;
        return;
    end
    
    % Check whether this model is suitable for fast sensitivities
    for m = 1 : numel( ar.model )
        arConservedPools(m);
        if numel( ar.model(m).pools.states ) > 0
            error( 'There are conserved moieties in model. Please run arReduce after arLoadModel and before loading data.' );
        end
    end
    
    % We're good!
    arFprintf(2, 'Enabled fast sensitivity computation\n');
    ar.config.turboSSSensi = 1;
end