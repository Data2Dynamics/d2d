% Attempts to find initial conditions such that the rhs for a specific
% condition equals the null vector.
%
%   [xnew, S] arFindRoots( jm, jc, condis )
%
% Usage:
%       jm      - Model number
%       jc      - Condition number
%       condis  - Condition (either "condition" or "ss_condition")
%
% Returns:
%       xnew    - Determined initial condition
%       S       - Sensitivities at determined point
%
% Note: This is an internal function

function [xnew, S] = arFindRoots(jm, jc, condis, useConserved)

    global ar;
    debug = 0;
    tolerance = .01 * ar.config.eq_tol;
    
    if nargin < 1
        jm = 1;
    end
    if nargin < 2
        jc = 1;
    end
    if nargin < 3
        condis = 'ss_condition';
    end
    if nargin < 4
        useConserved = 1;
    end
    
    % Determine a reference sensitivity for debugging purposes
    if ( debug )
        arCheckCache(1);
        arSimu(true, false, true);
        ar.model(jm).ss_condition(jc).sxFineSimu(end,:,:)
    end
    
    nS = length( ar.model(jm).x );
    
    % Grab initial x0 based on model parameters
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, 'ss_condition', 'ss_threads', 1);
    x0 = ar.model(jm).ss_condition(jc).xFineSimu(1,:);    
    
    if ( useConserved )
        % Get conserved pools (need independent initials eventually)
        dependent = null(ar.model(jm).N.', 'rational');
    
        % Set up the objective function for lsqnonlin
        fn = @(x)merit( x, jm, jc );
    end

    % Estimate initials in steady state
    opts            = optimset('TolFun', tolerance*tolerance, 'Display', 'Off' );
    [xnew, resnorm] = lsqnonlin( fn, x0, 0*x0, [], opts );
    
    % Calculate sensitivities via implicit function theorem
    dfdx = ar.model.N * ar.model(jm).(condis)(jc).dvdxNum;
    dfdp = ar.model.N * ar.model(jm).(condis)(jc).dvdpNum;
    S    = -pinv(dfdx)*dfdp;
    
    if ( resnorm > tolerance )
        warning( 'Failure to converge when rootfinding for model %d, condition %d', jm, jc );
    end
    
    % Remove the override after determination
    ar.model(jm).ss_condition(jc).x0_override = [];
end

% dxdts are squared to generate minimum for small dxdt
function res = merit(x0, jm, jc)
    global ar;
    ar.model(jm).ss_condition(jc).x0_override = x0;
    
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, 'ss_condition', 'ss_threads', 1);
    res = ar.model(jm).ss_condition(jc).dxdt.*ar.model(jm).ss_condition(jc).dxdt;
end