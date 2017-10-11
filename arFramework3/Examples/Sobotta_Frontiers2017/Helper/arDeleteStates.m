% This function can be used to remove states.
% WARNING: This function can only be used immediately after loading a model 
% and *before* loading any data. Using the function after loading data will
% result in model corruption and potentially segmentation faults.

function arDeleteStates( m, removedStates )

    global ar;

    notFound = setdiff( removedStates, ar.model(m).x );
    if ( length( notFound ) > 0 )
        fprintf( 'Did not find the following states: ' );
        fprintf( '%s ', notFound{:} );
    end

    % Remove states
    statesToRemove = ismember(ar.model(m).x, removedStates);
    ar.model(m).x(statesToRemove) = [];
    ar.model(m).fx(statesToRemove) = [];
    ar.model(m).qPlotX(statesToRemove) = [];
    ar.model(m).qPositiveX(statesToRemove) = [];
    ar.model(m).xNames(statesToRemove) = [];
    ar.model(m).xUnits(statesToRemove,:) = [];
    ar.model(m).N(statesToRemove, :) = [];
    ar.model(m).cLink(statesToRemove) = [];
    
    % Remove reactions that are no longer used
    redundantReactions = max(abs(ar.model(m).N), [], 1) < 1;
    ar.model(m).v(redundantReactions) = [];
    ar.model(m).fv(redundantReactions) = [];
    ar.model(m).qPlotV(redundantReactions) = [];
    if ( isfield( ar.model(m), 'vNames' ) )
        ar.model(m).vNames(redundantReactions) = [];
    end
    ar.model(m).vUnits(redundantReactions,:) = [];
    ar.model(m).N(:, redundantReactions) = [];
    ar.model(m).fv_ma_reverse_pbasename(redundantReactions) = [];
    ar.model(m).fv_source(redundantReactions) = [];
    ar.model(m).fv_sourceCoeffs(redundantReactions) = [];
    ar.model(m).fv_target(redundantReactions) = [];
    ar.model(m).fv_targetCoeffs(redundantReactions) = [];
    
    % Remove sources and targets which no longer exist
    for js = 1 : length( ar.model(m).fv_source )
        if ismember( ar.model(m).fv_source{js}, removedStates )
            ar.model(m).fv_source{js} = {};
            ar.model(m).fv_sourceCoeffs{js} = [];
        end
    end
    for js = 1 : length( ar.model(m).fv_target )
        if ismember( ar.model(m).fv_target{js}, removedStates )
            ar.model(m).fv_target{js} = {};
            ar.model(m).fv_targetCoeffs{js} = [];
        end
    end
    
    % remove states from rxn eqs
    symStates = sym( removedStates );
    for jv=1:length(ar.model(m).fv)
        exp = sym( ar.model(m).fv{jv} );
        exp = subs( exp, symStates, zeros(size(symStates)) );
        ar.model(m).fv{jv} = char(exp);
    end
    
    % remove states from rxn eqs
    for jx=1:length(ar.model(m).fx)
        exp = sym( ar.model(m).fx{jx} );
        exp = subs( exp, symStates, zeros(size(symStates)) );
        ar.model(m).fx{jx} = char(exp);
    end
    
    % remove derived quantities that depend on the removed states
    removeDerived = false( 1, numel( ar.model(m).fz ) );
    for jd=1:length(ar.model(m).fz)
        vars = symvar(ar.model(m).fz{jd});
        if max(ismember(vars, removedStates))
            removeDerived(jd) = true;
        end
    end
    removedDerived = ar.model(m).z(removeDerived);
    ar.model(m).z(removeDerived) = [];
    ar.model(m).qPlotZ(removeDerived) = [];
    
    % remove observables that depend on the removed states
    removeObs = false( 1, numel( ar.model(m).fy ) );
    removedQuantities = union( removedStates, removedDerived );
    for jo=1:length(ar.model(m).fy)
        vars = symvar(ar.model(m).fy{jo});
        if max(ismember(vars, removedQuantities))
            removeObs(jo) = true;
        end
    end
    ar.model(m).y(removeObs) = [];
    ar.model(m).fy(removeObs) = [];
    ar.model(m).fystd(removeObs) = [];
    ar.model(m).yUnits(removeObs,:) = [];
    ar.model(m).yNames(removeObs) = [];
    ar.model(m).logfitting(removeObs) = [];
    ar.model(m).logplotting(removeObs) = [];
    ar.model(m).normalize(removeObs) = [];

    % reconstruct dynamic parameters
    ar.model(m).px0 = strcat('init_', ar.model.x);
    new_px = ar.model.pc;
    
    % inputs
    varlist = cellfun(@symvar, ar.model(m).fu, 'UniformOutput', false);
    ar.model(m).pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible
        
    % parameters from the flux expressions
    ar.model(m).pvs = cell(size(ar.model(m).fv));
    ar.model(m).pv = {};
    for jv=1:length(ar.model(m).fv)
        varlist = symvar(ar.model(m).fv{jv});
        ar.model(m).pvs{jv} = setdiff(varlist, union(ar.model(m).t, union(ar.model(m).x, ar.model(m).u))); %R2013a compatible
        ar.model(m).pv = union(ar.model(m).pv, ar.model(m).pvs{jv});
    end
    
    % derived variables parameters
    varlist = cellfun(@symvar, ar.model(m).fz, 'UniformOutput', false);
    ar.model(m).pz = setdiff(setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}), union(ar.model(m).x, union(ar.model(m).u, ar.model(m).z))); %R2013a compatible
    ar.model(m).px = union( union( union( ar.model(m).pv, new_px), ar.model(m).px0 ), ar.model(m).pz ); %R2013a compatible
        
    % observation parameters
    varlist = cellfun(@symvar, ar.model(m).fy, 'UniformOutput', false);
    ar.model(m).py = setdiff(setdiff(vertcat(varlist{:}), union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ar.model(m).z)), {ar.model(m).t, ''});
    varlist = cellfun(@symvar, ar.model(m).fystd, 'UniformOutput', false);
    ar.model(m).pystd = setdiff(vertcat(varlist{:}), union(union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ar.model(m).z), ar.model(m).y), ar.model(m).t));
    
    % Assemble new parameter vector
    new_p = union(union(union(ar.model(m).px, ar.model(m).pu), ar.model(m).py), ar.model(m).pystd);        
    
    % Remove parameters which disappeared from the model
    discardedParameters = ~ismember(ar.model(m).p, new_p);
    ar.model(m).fp(discardedParameters) = [];
    ar.model(m).p(discardedParameters) = [];
        
    % Reconstruct condition variables
    varlist = cellfun(@symvar, ar.model(m).fp, 'UniformOutput', false);
    ar.model(m).pcond = setdiff(setdiff(setdiff(vertcat(varlist{:}), ar.model(m).p), ar.model(m).x), ar.model(m).u); %R2013a compatible
    
end