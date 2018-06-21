function arConservedPools(jm, showPools)
    
    global ar;
    
    if ( nargin < 2 )
        showPools = 0;
    end
    nx = length( ar.model(jm).x );
    
    % Try treating the compartment sizes as numbers, otherwise fall back to
    % symbolic math.
    try
        pc = cellfun(@str2num, ar.model.pc);
        wantSyms = 0;
        always = @(x)x;
    catch
        pc = sym(ar.model.pc);
        assume( pc > 0 );
        always = @(x)condAlways(x);
        wantSyms = 1;
    end
    
    pc = pc(ar.model(jm).cLink);
    N  = ar.model(jm).N;
    
    dependent = null( N.', 'r' );
    [~,i]=sort(sum(dependent),'descend');
    dependent = dependent(:, i);
    
    IDs = repmat(1:size(dependent,1), size(dependent,2), 1).';
    groups = zeros( size(IDs) );
    if ( wantSyms )
        groups = sym(groups);
    end

    % Account for compartment sizes
    % selected = dependent .* repmat(pc.', 1, size(dependent,2));
    dependent = dependent .* repmat(pc.', 1, size(dependent,2));
    selected = dependent;
    
    states = zeros( 1, size( dependent, 2 ) );
    for js = 1 : size( dependent, 2 )
        removedState = find( always(abs(selected) > 0), 1 );
        curCoeff = selected(removedState, 1);
        states(js) = IDs(removedState,1);
        IDs(removedState, :) = [];
        selected(removedState, :) = [];
        
        % List of other species involved in this relation
        otherSpecies = IDs(:,1);
        otherSpecies = otherSpecies( always( abs( selected(:,1) ) > 0 ) );
        stoich       = selected( always( abs( selected(:,1) ) > 0 ) );
        
        % Store coefficients of other species involved in the conserved pool
        groups( otherSpecies, js ) = - stoich / curCoeff;
        
        IDs(:, 1) = [];
        selected(:, 1) = [];
    end
    
    % Mapping from reduced model states to full states
    mapping = eye(nx);
    if ( wantSyms )
        mapping = sym( mapping );
    end
    for js = 1 : length( states )
        mapping(states(js),:) = groups(:,js);
    end
    mapping(:, states) = [];
    
    % Matrix to compute totals for respective species
    totalMapping = zeros(nx);
    if ( wantSyms )
        totalMapping = sym( totalMapping );
    end
    totalMapping(states, :) = dependent.';
    
    ar.model(jm).pools.dependent  = dependent;
    ar.model(jm).pools.others     = groups;
    ar.model(jm).pools.states     = states;
    ar.model(jm).pools.usedStates = setdiff( 1:length(ar.model(jm).x), states );
    ar.model(jm).pools.mapping    = mapping;
    ar.model(jm).pools.totalMap   = totalMapping;
    ar.model(jm).pools.isSymbolic = wantSyms;
    
    if ( showPools )
        disp( 'Reparameterized model:' );
        for jg = 1 : size( groups, 2 )
            c = '';
            for js = 1 : size( groups, 1 )
                if ( abs( groups(js, jg) ) > 0 )
                    c = sprintf( '%s + %d %s ', c, groups(js, jg), ar.model(jm).x{js} );
                end
            end
            fprintf( '%s = total_%d %s\n', ar.model(jm).x{states(jg)}, jg, c );
        end
    end 
end

function x = condAlways(x)
    if ~islogical(x)
        x = isAlways(x);
    end
end