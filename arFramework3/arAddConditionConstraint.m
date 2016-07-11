% function arAddConditionConstraint( m, c1, c2, t, sd, states, relative )
%   
%   Mandatory arguments:
%       m           Model within which to constrain
%       c1          Condition 1
%       c2          Condition 2
%       t           Time points to use in the constraint
%   Optional arguments:
%       states      States to constrain (default = 1:nStates)
%       relative    Constrain on log10-scale? (default = 1)
%
% Warning: This function invokes arLink, which may lead you to lose
% simulation results.
%
% Warning: This functionality is still being developed and the interface is 
% subject to change. Please be aware that backward compatibility for this
% function is not ensured at this time.

function arAddConditionConstraint( m, c1, c2, t, sd, states, relative )

    global ar;
    
    strct.m         = m;
    strct.c1        = c1;
    strct.c2        = c2;
    strct.t         = t;
    strct.sd        = sd;
    
    if ( ~exist( 'states', 'var' ) )
        strct.states = 1:length(ar.model(m).x);
    else
        strct.states = states;
    end
    
    if ( ~exist( 'relative', 'var' ) )
        strct.relative = 1;
    else
        strct.relative = relative;
    end

    if ~isfield( ar, 'conditionconstraints' )
        ar.conditionconstraints = strct;
    else
        ar.conditionconstraints(end+1) = strct;
    end

    % Add necessary time points to the conditions
    arAddEvent(m, c1, t );
    arAddEvent(m, c2, t );
    arLink;
    
end