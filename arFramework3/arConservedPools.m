function arConservedPools(jm, showPools)
    
    global ar;
    
    if ( nargin < 2 )
        showPools = 0;
    end
    
    N = ar.model(jm).N;
    dependent = null( N.', 'r' );
    
    IDs = repmat(1:size(dependent,1).', size(dependent,2), 1).';
    groups = zeros( size(IDs) );
    
    selected = dependent;
    for js = 1 : size( dependent, 2 )
        removedState = find( selected == 1, 1 );
        curCoeff = selected(removedState, 1);
        states(js) = IDs(removedState,1);
        IDs(removedState, :) = [];
        selected(removedState, :) = [];
        
        % List of other species involved in this relation
        otherSpecies = IDs(:,1);
        otherSpecies = otherSpecies(  abs( selected(:,1) ) > 0 );
        stoich       = selected(  abs( selected(:,1) ) > 0 );
        
        % Store coefficients of other species involved in the conserved pool
        groups( otherSpecies, js ) = - stoich / curCoeff;
        
        IDs(:, 1) = [];
        selected(:, 1) = [];
    end
    
    % Mapping from reduced model states to full states
    mapping = eye(length(ar.model.x));    
    for js = 1 : length( states )
        mapping(states(js),:) = groups(:,js);
    end
    mapping(:, states) = [];
    
    % Matrix to compute totals for respective species
    totalMapping = zeros(length(ar.model.x));
    totalMapping(states, :) = dependent.';
    
    ar.model(jm).pools.dependent = dependent;
    ar.model(jm).pools.others    = groups;
    ar.model(jm).pools.states    = states;
    ar.model(jm).pools.mapping   = mapping;
    ar.model(jm).pools.totalMap  = totalMapping;
    
    if ( showPools )
        disp( 'Reparameterized model:' );
        for jg = 1 : size( groups, 2 )
            c = '';
            for js = 1 : size( groups, 1 )
                if ( abs( groups(js, jg) ) > 0 )
                    c = sprintf( '%s + %d %s ', c, groups(js, jg), ar.model.x{js} );
                end
            end
            fprintf( '%s = total_%d %s\n', ar.model.x{states(jg)}, jg, c );
        end
    end