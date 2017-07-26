% This function needs to be called before simulating with the sensitivity
% subset option activated. It prepares some data structures required for
% determining which sensitivities to simulate.

function arSubsetSensi()

    global ar;

    % Establish the correct sensitivity mappings
    for jm = 1 : numel( ar.model )
        for jc = 1 : numel( ar.model(jm).condition )
            ar.model(jm).condition(jc).sensIndices = find( ar.qFit( ar.model(jm).condition(jc).pLink ) );
            ar.model(jm).condition(jc).backwardIndices = zeros( sum( ar.model(jm).condition(jc).pLink ), 1 );
            ar.model(jm).condition(jc).backwardIndices( ar.model(jm).condition(jc).sensIndices ) = 1 : numel( ar.model(jm).condition(jc).sensIndices );
            
            % Convert to usable datatype
            ar.model(jm).condition(jc).sensIndices = int32( ar.model(jm).condition(jc).sensIndices - 1 );
            ar.model(jm).condition(jc).backwardIndices = int32( ar.model(jm).condition(jc).backwardIndices - 1 );
        end
    end
    
    % Establish the correct sensitivity mappings
    for jm = 1 : numel( ar.model )
        if ( isfield( ar.model(jm), 'ss_condition' ) )
            for jc = 1 : numel( ar.model(jm).ss_condition )
                ar.model(jm).ss_condition(jc).sensIndices = find( ar.qFit( ar.model(jm).ss_condition(jc).pLink ) );
                ar.model(jm).ss_condition(jc).backwardIndices = zeros( sum( ar.model(jm).ss_condition(jc).pLink ), 1 );
                ar.model(jm).ss_condition(jc).backwardIndices( ar.model(jm).ss_condition(jc).sensIndices ) = 1 : numel( ar.model(jm).ss_condition(jc).sensIndices );
    
                % Convert to usable datatype
                ar.model(jm).ss_condition(jc).sensIndices = int32( ar.model(jm).ss_condition(jc).sensIndices - 1 );
                ar.model(jm).ss_condition(jc).backwardIndices = int32( ar.model(jm).ss_condition(jc).backwardIndices - 1 );
            end
        end
    end
        