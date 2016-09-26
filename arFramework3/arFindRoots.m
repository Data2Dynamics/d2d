% Work in progress... do not use yet!

function arFindRoots()

    global ar;
    debug = 0;
    jm = 1;
    jc = 1;
    tolerance = .01 * ar.config.eq_tol;
    
    condis = 'ss_condition';
    
    %%% Get a reference sensitivity for now
    if ( debug )
        arCheckCache(1);
        arSimu(true, false, true);
        ar.model(jm).ss_condition(jc).sxFineSimu(end,:,:)
    end
    
    nS = length( ar.model(jm).x );
    
    % Get conserved pools (need independent initials eventually)
    dependent = null(ar.model(jm).N.', 'rational');
    
    % Do not simulate, but grab initial condition that is set in the model
    ar.model(jm).(condis)(jc).dvdxNum = 0 * ar.model(jm).(condis)(jc).dvdxNum;
    ar.model(jm).(condis)(jc).dvduNum = 0 * ar.model(jm).(condis)(jc).dvduNum;
    ar.model(jm).(condis)(jc).dvdpNum = 0 * ar.model(jm).(condis)(jc).dvdpNum;
    
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, 'ss_condition', 'ss_threads', 1);
    x0 = ar.model(jm).ss_condition(jc).xFineSimu(1,:);
    
    % Set up the objective function for lsqnonlin
    fn = @(x)merit( x, jm, jc );

    % Estimate initials in steady state
    opts            = optimset('TolFun', tolerance*tolerance, 'Display', 'Off' );
    [xnew resnorm]  = lsqnonlin( fn, x0, 0*x0, [], opts );
    
    % Calculate sensitivities via implicit function theorem
    dfdx = ar.model.N * ar.model(jm).(condis)(jc).dvdxNum;
    dfdu = ar.model.N * ar.model(jm).(condis)(jc).dvduNum;
    dfdp = ar.model.N * ar.model(jm).(condis)(jc).dvdpNum;
    S = -inv(dfdx)*dfdp;
    
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