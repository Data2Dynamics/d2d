% This function checks whether any of the integration settings changed. It
% is used in arSimu to determine whether sensitivity equations should be
% simulated or not.
%
% function upToDate = arCheckCache( invalidate )
% Set invalidate to 1 to invalidate the current cache. This forces
% resimulation of the sensitivities on the next simulation.

function invalidate = arCheckCache( invalidate )
    global ar;

    fields = {  'useParallel', 'useJacobian', 'useSparseJac', 'useSensiRHS', ...
                'atolV', 'atolV_Sens', 'rtol', 'atol', 'maxsteps', 'maxstepsize', ...
                'fiterrors', 'ssa_min_tau', 'ssa_runs', 'useMS', 'useEvents', ...
                'max_eq_steps', 'init_eq_step', 'eq_step_factor', 'eq_tol', 'sensiSkip' };

    if ( nargin < 1 )
        invalidate = 0;
    end
    
    % Check if fields are available in ar.config
    fields = checkArConfigFields( fields );

    % Check whether the config settings were the same for the previous
    % simulation
    if ( ~invalidate )
        invalidate = ~checkCacheConfigFields( fields );
    end
    
    % Invalid cache, remove both cached fine and exp field
    if ( invalidate )
        ar.cache.fine   = nan(size(ar.p));
        ar.cache.exp    = nan(size(ar.p));
        
        % Set the cache to the current cache
        setCacheConfigFields( fields );
    end
    
end

function fields = checkArConfigFields( fields )
    global ar;
    
    remove = zeros(1,length( fields ));
    for a = 1 : length( fields )
        if ~isfield( ar.config, fields{a} )
            remove(a) = 1;
        end
    end
    fields = fields(~remove);
end

function setCacheConfigFields( fields )
    global ar;

    for a = 1 : length( fields )
        ar.cache.(fields{a}) = ar.config.(fields{a});
    end
end

function valid = checkCacheConfigFields( fields )
    global ar;

    if ~isfield( ar, 'cache' )
        valid = 0;
        return;
    end
    
    for a = 1 : length( fields )
        if ~isequal( ar.cache.(fields{a}), ar.config.(fields{a}) );
            ar.cache.valid = 0;
            valid = 0;
            return
        end
    end
   
    valid = 1;
    return;
end