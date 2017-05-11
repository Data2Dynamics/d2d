% This function checks whether any of the integration settings changed. It
% is used in arSimu to determine whether sensitivity equations should be
% simulated or not.
%
% function upToDate = arCheckCache( invalidate )
% Set invalidate to 1 to invalidate the current cache. This forces
% resimulation of the sensitivities on the next simulation.

function invalidate = arCheckCache( invalidate )
    global ar;

    if ( nargin < 1 )
        invalidate = 0;
    end

    if ( (~isfield( ar, 'info' )) || (~isfield( ar.info, 'arFormatVersion' )) || ( ar.info.arFormatVersion ~= arInitFields() ) )
        arFprintf( 0, 'Adding missing fields (updating ar struct from older ar version)\n' );
        ar = arInitFields(ar);
    end
    
    % Always resimulate
    if ( ar.config.useCache == 0 )
        invalidate = 1;
    end
    
    fields = fieldnames( ar.config );

    % Check whether the config settings were the same for the previous
    % simulation
    if ( ~invalidate )
        invalidate = ~checkCacheConfigFields( fields );

        % For backward compatibility
        if ( isfield( ar, 'cache' ) )
            if ( ~isfield( ar.cache, 'expSensi' ) )
                invalidate = 1;
            end
        else
            invalidate = 1;
        end
    end
    
    % Invalid cache, remove both cached fine and exp field
    if ( invalidate )
        ar.cache.fine           = nan(size(ar.p));
        ar.cache.exp            = nan(size(ar.p));
        ar.cache.fineSensi      = nan;
        ar.cache.expSensi       = nan;
               
        % Set the cache to the current cache
        setCacheConfigFields( fields );
    end
    
end

function setCacheConfigFields( fields )
    global ar;

    for a = 1 : length( fields )
        if ( isnumeric( ar.config.(fields{a}) ) )
            ar.cache.(fields{a}) = ar.config.(fields{a}) + 0;
        else
            ar.cache.(fields{a}) = ar.config.(fields{a});
        end
    end
end

function valid = checkCacheConfigFields( fields )
    global ar;

    if ~isfield( ar, 'cache' )
        valid = 0;
        return;
    end
    
    valid = 1;
    for a = 1 : length( fields )
        try
            if ~isequal( ar.cache.(fields{a}), ar.config.(fields{a}) )
                if ( ~( isnan( ar.cache.(fields{a}) ) == isnan( ar.config.(fields{a}) ) ) )
                    ar.cache.valid = 0;
                    valid = 0;
                    return;
                end
            end
        catch
            valid = 0;
            return;
        end
    end
end