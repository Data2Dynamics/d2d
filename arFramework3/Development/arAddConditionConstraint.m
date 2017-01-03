% function arAddConditionConstraint( m, c1, c2, t, sd, states, relative )
%   
%   Mandatory arguments:
%       m1          Model 1, within which to constrain condition 1
%       c1          Condition 1
%       m2          Model 2, within which to constrain condition 2
%       c2          Condition 2
%       t           Time points to use in the constraint
%   Optional arguments:
%       states      Cell array of state names to constrain (default = 'all'  
%                   which matches all matching states)
%       relative    Constrain on log10-scale? (default = 1)
%
% Warning: This function invokes arLink, which may lead you to lose
% simulation results.
%
% Warning: This functionality is still being developed and the interface is 
% subject to change. Please be aware that backward compatibility for this
% function is not ensured at this time.

function arAddConditionConstraint( m1, c1, m2, c2, t, sd, states, relative, silent )

    global ar;
    
    if ( ~exist( 'silent', 'var' ) )
        silent = 0;
    end
    
    strct.m1        = m1;
    strct.c1        = c1;
    strct.m2        = m2;
    strct.c2        = c2;
    strct.t         = t;
    strct.sd        = sd;
    
    if ( ~exist( 'states', 'var' ) || strcmp( states, 'all' ) )
        states = union( ar.model(m1).x, ar.model(m2).x );
    end
    
    % Find which states exist in both models
    sharedStates = intersect( intersect( ar.model(m1).x, states ), intersect( ar.model(m2).x, states ) );
    
    % Find IDs for linkage
    strct.states1 = find( ismember( ar.model(m1).x, sharedStates ) );
    strct.states2 = find( ismember( ar.model(m2).x, sharedStates ) );
    
    strlen = max( [ cellfun( @length, ar.model(m1).x ), cellfun( @length, ar.model(m2).x ) ] );
    
    if ( ~silent )
        fprintf( 'Constrained states:\n' );
        for a = 1 : length( strct.states1 )
            fprintf( 'm%d.%s [%d] <=> m%d.%s [%d]\n', m1, arExtendStr( ar.model(m1).x{strct.states1(a)}, strlen ), strct.states1(a), ...
                        m2, arExtendStr( ar.model(m2).x{strct.states2(a)}, strlen ), strct.states1(a) );
        end
    end
    
    if ( ~exist( 'relative', 'var' ) )
        strct.relative = 1;
    else
        strct.relative = relative;
    end
    
    strct.tLink1    = [];
    strct.tLink2    = [];
    
    if ~isfield( ar, 'conditionconstraints' )
        ar.conditionconstraints = strct;
    else
        ar.conditionconstraints(end+1) = strct;
    end

    % Add necessary time points to the conditions
    arAddEvent(m1, c1, t );
    arAddEvent(m2, c2, t );
    arLink;
    
end