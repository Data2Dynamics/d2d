% function arReduce( [m], [keepTotals] )
%
% This function automatically reduces the model by removing one state for
% every conserved moiety. It needs to be called right after arLoadModel if 
% this method is to be used. Before loading data and before compilation!
%
%   m               Model index                     [1]
%   keepTotals      Express model in totals?        [0]
%
% Example(s):
%   arInit;
%   arLoadModel(mymodel);
%   arReduce;
%   arLoadData etc.
%   arCompileAll;
%
% Note: arReduce will not change the model parameters. Internally, it
% reformulates to a formulation using 1 total pool and N-1 states per
% conserved moieties; but the model parameters will remain the same. The
% mapping from independent initials to totals is done internally (i.e. 
% by default you will *not* end up with a formulation in terms of totals).
%
% States which are removed are placed in the derived variables (ar.model.z),
% such that they can still be used in observation functions.
%
% This functionality can be used to greatly speed up models that use 
% pre-equilibrated steady states. In this case, simulation of the
% sensitivity equations can be avoided when equilibrating the model.
% This rapid equilibration can be toggled via the command arFastSensis.
%
% See also: arSteadyState, arFastSensis
%
function arReduce( m, keepTotals )
    global ar;
    
    if ( nargin < 1 )
        m = 1;
    end
    if ( nargin < 2 )
        keepTotals = 0;
    end
    
    if isfield( ar.model(m), 'condition' )
        error( 'This function should be called before loading any data' );
    end
    if isfield( ar.model(m), 'data' )
        error( 'This function should be called before loading any data' );
    end    
    if isfield( ar.model(m), 'reducedForm' )
        error( 'This model has already been reduced' );
    end
    ar.model(m).reducedForm = 1;
    
    % Fetch the conserved pools
    arConservedPools(m);
    
    % Are we dealing with symbolic compartments?
    if ( ar.model(m).pools.isSymbolic )
        always = @(x)condAlways(x);     
    else
        always = @(x)x;
    end
       
    arFprintf( 2, 'Pools available for reduction:\n' );
    for i = 1 : size(ar.model(m).pools.dependent, 2)
        arFprintf( 2, '%d: %s\n', i, sprintf( '%s ', ar.model(m).x{always(abs(ar.model(m).pools.dependent(:,i))>0)} ) );
    end
    
    if ( isempty( ar.model(m).pools.states ) )
        warning( 'Nothing to reduce (No conserved moieties / dfdx already full rank)' );
        ar.model(m).reducedForm = 0;
    end
    
    while( ~isempty( ar.model(m).pools.states ) )
        % Add replacement rule for where the state appears in equations
        [removedState, removalStruct] = genReplacementRule( m, 1, always, keepTotals );
        if isfield( ar.model(m), 'removedStates' )
            ar.model(m).removedStates(end+1) = removalStruct;
        else
            ar.model(m).removedStates = removalStruct;
        end
        ar.model(m).px{end+1}   = removalStruct.totalVariable;
        ar.model(m).p{end+1}    = removalStruct.totalVariable;
        ar.model(m).fp{end+1}   = removalStruct.totalVariable;
        
        % Since unions are typically used in the code to merge the parameters
        % and this sorts them, we need to make sure we keep our parameters
        % sorted
        [ar.model(m).p, I]      = sort( ar.model(m).p );
        ar.model(m).fp          = ar.model(m).fp(I);
        
        ar.model(m).z{end+1}    = removalStruct.p;
        ar.model(m).fz{end+1}   = removalStruct.fp;
        
        % Substitute state variable in rate equations
        fvtemp = arSubs( arSym( ar.model(m).fv ), removalStruct.sym.p, removalStruct.sym.fp );
        fv = cell( 1, numel( fvtemp ) );
        for a = 1 : numel( fvtemp )
            fv{a} = char( fvtemp(a) );
        end
        ar.model(m).fv = fv.';
        
        % Remove the states from the stoichiometry
        ar.model(m).N( removedState, : )       = [];
        ar.model(m).Cm( removedState, : )      = [];
        ar.model(m).Cm_par( removedState, : )  = [];

        name = ar.model(m).x( removedState );
        for i = 1 : numel( ar.model(m).fv_source )
            idx = find( strcmp( ar.model(m).fv_source{i}, name ) );
            if ~isempty( idx )
                ar.model(m).fv_source{i}(idx) = [];
                ar.model(m).fv_sourceCoeffs{i}(idx)=[];
            end
        end
        for i = 1 : numel( ar.model(m).fv_target )
            idx = find( strcmp( ar.model(m).fv_target{i}, name ) );
            if ~isempty( idx )
                ar.model(m).fv_target{i}(idx) = [];
                ar.model(m).fv_targetCoeffs{i}(idx)=[];
            end
        end        
        
        % Remove the states from the initial conditions
        ar.model(m).px0( removedState )        = [];
        ar.model(m).zUnits( end + 1, : )       = ar.model(m).xUnits( removedState, : );
        if isfield( ar.model(m), 'zNames' )
            ar.model(m).zNames( end + 1, : )       = ar.model(m).xNames( removedState );
        end
        ar.model(m).xUnits( removedState, : )  = [];
        ar.model(m).cLink( removedState )      = [];
        ar.model(m).x( removedState )          = [];
        ar.model(m).xNames( removedState )     = [];
        ar.model(m).qPlotX( removedState )     = [];
        ar.model(m).qPlotZ( end + 1 )          = 1;        
        
        % Update conserved pools one at a time
        arConservedPools(m);
    end
    
    % We've completed the iterative reductions, we can now look for the new
    % states in the new model. This is not necessary, but useful if we are
    % to diagnose problems at a later stage.
    %findNewStates(m);
      
    % Regenerate equations
    tmpfx = (arSym(ar.model(m).N).*arSym(ar.model(m).Cm)) * arSym(ar.model(m).fv);
    tmpfx_par = (arSym(ar.model(m).N).*arSym(ar.model(m).Cm_par)) * arSym(ar.model(m).fv);

    for j=1:length(ar.model(m).x) % for every species j
        if ~isempty(tmpfx)
            ar.model(m).fx{j}       = char(tmpfx(j));
            ar.model(m).fx_par{j}   = char(tmpfx_par(j));
        else
            ar.model(m).fx{j}       = char('0');
            ar.model(m).fx_par{j}   = char('0');
        end
    end
    
    % If nothing was done, make sure this field exists but is empty
    % We did still check whether the dfdx was full rank, so it should be OK
    % to use implicit methods.
    if ~isfield( ar.model(m), 'removedStates' )
        ar.model(m).removedStates = [];
    end
end

% Generate replacement rules
function [removedState, removalStruct] = genReplacementRule( m, selected, always, keepTotals )
    global ar;
    
    removedState = ar.model(m).pools.states(selected);
    p = ar.model(m).x{ removedState };
    
    % Make unique variable which holds the total pool size
    totalVariable = [ '__total_', ar.model(m).x{always(abs(ar.model(m).pools.dependent(:,selected))>0)} ];
    otherStates = {};
    
    % The total has to be divided by the stoichiometric contribution to the
    % total pool of the state we are removing (the conservation relations are 
    % stored in ar.model(m).pools.dependent)
    if ( ~isnumeric( ar.model(m).pools.dependent(removedState, selected) ) )
        fp = sprintf( '(%s) * %s', char(1 / ar.model(m).pools.dependent(removedState, selected)), totalVariable );
    else
        fp = sprintf( '%g * %s', 1 / ar.model(m).pools.dependent(removedState, selected), totalVariable );
    end
    
    % Grab the other variables involved in this pool. These are stored in 
    % ar.model(m).pools.others; along with their respective stoichiometry.
    other = find( abs( ar.model(m).pools.others( :, selected ) ) > 0 );
    for js = 1 : numel( other )
        if ( ~isnumeric( ar.model(m).pools.others(other(js), selected) ) )
            fp = sprintf( '%s + %s * %s', fp, char( ar.model(m).pools.others(other(js), selected) ), ar.model(m).x{other(js)} );
        else
            fp = sprintf( '%s + %g * %s', fp, ar.model(m).pools.others(other(js), selected), ar.model(m).x{other(js)}  );
        end
        
        otherStates{end+1} = { ar.model(m).pools.others(other(js)), ar.model(m).x{other(js)} }; %#ok
    end
    
    % Determine which initial conditions will contribute to this total pool
    % such that we can compute the value for this parameter in each of the
    % conditions later
    poolIDs = find(always(abs(ar.model(m).pools.dependent(:,selected))>0));
    totalPoolStates = '';
    for js = 1 : numel( poolIDs )
        if ( ~isnumeric( ar.model(m).pools.dependent(poolIDs(js), selected) ) )
            totalPoolStates = sprintf( '%s + %s * init_%s', totalPoolStates, char( ar.model(m).pools.dependent(poolIDs(js), selected) ), ar.model(m).x{poolIDs(js)}  );
        else
            totalPoolStates = sprintf( '%s + %g * init_%s', totalPoolStates, ar.model(m).pools.dependent(poolIDs(js), selected), ar.model(m).x{poolIDs(js)}  );
        end
    end
    
    % Store the fetched and generated results
    removalStruct.expressInTotal = keepTotals;
    removalStruct.p = p;
    removalStruct.fp = fp;
    removalStruct.initial = ['init_' p];
    removalStruct.totalVariable = totalVariable;
    removalStruct.totalPoolStates = totalPoolStates;
    removalStruct.otherStates = otherStates;
    removalStruct.sym.p = arSym(p);
    removalStruct.sym.fp = arSym(fp);
    removalStruct.sym.totalVariable = arSym(totalVariable);
    removalStruct.sym.totalPoolStates = arSym(totalPoolStates);
    removalStruct.sym.initial = arSym(removalStruct.initial);
end

% Convert the other states from their name to their indices as preparation
% for the compilation process
function findNewStates(m)
    global ar;
    for jrs = 1 : numel( ar.model(m).removedStates )
        for jrc = 1 : numel( ar.model(m).removedStates(jrs).otherStates )
            coefficient = ar.model(m).removedStates(jrs).otherStates{jrc}{1};
            newStateID = find( strcmp( ar.model(m).x, ar.model(m).removedStates(jrs).otherStates{jrc}{2} ) );
            
            if isfield( ar.model(m).removedStates(jrs), 'removedStateCoeffs' )
                ar.model(m).removedStates(jrs).removedStateCoeffs(end+1) = coefficient;
                ar.model(m).removedStates(jrs).stateIDs(end+1) = newStateID;
            else
                ar.model(m).removedStates(jrs).removedStateCoeffs = coefficient;
                ar.model(m).removedStates(jrs).stateIDs = newStateID;
            end
        end
    end
end

function x = condAlways(x)
    if ~islogical(x)
        x = isAlways(x);
    end
end