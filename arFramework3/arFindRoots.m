% Work in progress... do not use yet!

function arFindRoots()

    jm = 1;
    jc = 2;
    global ar;
    
    %%% Get a reference sensitivity for now
    arCheckCache(1);
    arSimu(true, true, true);
    ar.model(jm).ss_condition(jc).sxFineSimu(end,:,:)
    %%%
    
    nS = length( ar.model(jm).x );
    
    % Get conserved pools (need independent initials eventually)
    dependent = null(ar.model(jm).N.', 'rational');
    
    % Do not simulate, but grab initial condition that is set in the model
    feval(ar.fkt, ar, true, false, true, false, 'ss_condition', 'ss_threads', 1);
    x0 = ar.model(jm).ss_condition(jc).xFineSimu(1,:);
    
    % Set up the objective function for lsqnonlin
    fn = @(x)merit( x, jm, jc );

    % Estimate initials in steady state
    opts = optimset('TolFun', 1e-8*1e-8, 'Display', 'Off' );
    xnew = lsqnonlin( fn, x0, 0*x0, [], opts );
    
    
    %feval(ar.fkt, ar, true, ar.config.useSensis && sensi, dynamics, false, 'ss_condition', 'ss_threads');
        
    % Calculate sensitivities via implicit function theorem (still something wrong)
    dfdx = ar.model.N * ar.model(jm).condition(jc).dvdxNum;
    dfdp = ar.model.N * ar.model(jm).condition(jc).dvdpNum;
    S = -inv(dfdx)*dfdp

    % Remove the override after determination
    ar.model(jm).ss_condition(jc).x0_override = [];
end

% dxdts are squared to generate minimum for small dxdt
function res = merit(x0, jm, jc)
    global ar;
    ar.model(jm).ss_condition(jc).x0_override = x0;
    
    feval(ar.fkt, ar, true, false, true, false, 'ss_condition', 'ss_threads', 1);
    res = ar.model(jm).ss_condition(jc).dxdt.*ar.model(jm).ss_condition(jc).dxdt;
end