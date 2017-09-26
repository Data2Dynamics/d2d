% Attempts to find initial conditions such that the rhs for a specific
% condition equals the null vector.
%
%   [xnew, S] arFindRoots( jm, jc, condis )
%
% Usage:
%       jm              - Model number
%       jc              - Condition number
%       condis          - Condition (either "condition" or "ss_condition")
%       useConserved    - Conserve conserved quantities
%       debug           - Provide debugging information
%       x0i             - Provide custom initial guess for x0
%
% Returns:
%       xnew            - Determined initial condition
%       S               - Sensitivities at determined point
%
% Note: This is an internal function

function [xnew, S, failedCheck] = arFindRoots(jm, jc, condis, useConserved, debug, x0i, newtonRaphson)

    global ar;
    tolerance = ar.config.eq_tol;
    
    if nargin < 1
        jm = 1;
    end
    if nargin < 2
        jc = 1;
    end
    if nargin < 3
        condis = 'ss_condition';
        threads = 'ss_threads';
    else
        if strcmp( condis, 'ss_condition' )
            threads = 'ss_threads';
        elseif strcmp( condis, 'condition' )
            threads = 'threads';
        else
            error( 'Condis has to be either condition or ss_condition' );
        end        
    end  
    if nargin < 4
        useConserved = 1;
    end
    if nargin < 5
        debug = 0;
        failedCheck = 0;
    end
    if nargin < 7
        newtonRaphson = 1;
    end

    % No need to incur the loss in speed when the model has been adequately
    % reduced
    if ( ar.model(jm).reducedForm )
        useConserved = 0;
    end
       
    % Determine a reference sensitivity for debugging purposes
    if ( debug )
        arCheckCache(1);
        arSimu(true, true, true);
        Sref = squeeze(ar.model(jm).(condis)(jc).sxFineSimu(end,:,:));
        xref = ar.model(jm).(condis)(jc).xFineSimu(end,:);
    end
    
    % Grab initial x0 based on model parameters
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, condis, threads, 1);
    x0 = ar.model(jm).(condis)(jc).xFineSimu(1,:);    
    
    if ( useConserved )
        warning('Custom:moieties','Model has not been reduced to its minimal form. Rootfinding with conserved moieties currently produces incorrect steady states when compartments are present. Please consider reducing the model first with arReduce (type help arReduce)')
        warning('off', 'Custom:moieties');
        
        % Get conserved pools (need independent initials eventually)
        if ( ~isfield( ar.model(jm), 'pools' ) )
            arConservedPools(jm);
        end
        
        % Compute total pools
        totals = ar.model(jm).pools.totalMap * x0.';
        
        % Get mapping from states with total pools substituted
        map = ar.model(jm).pools.mapping;
        
        % Compute reduced state vector
        x0( ar.model(jm).pools.states ) = [];
                
        % Set up the objective function for lsqnonlin
        fn = @(x)meritConserved( x, jm, jc, map, totals, condis, threads );
    else
        % Set up the objective function for lsqnonlin
        fn = @(x)merit( x, jm, jc, condis, threads );
    end
    
    if ( nargin < 6 ) || ( isempty( x0i ) )
        x0i = x0;
    end
    
    % Use Newton Raphson or lsqnonlin (former is typically much faster)? 
    if newtonRaphson
        if ( useConserved )
            % Have we not removed our conserved moieties, then go the slow route
            [xnew, maxDiff] = NewtonRaphson( fn, x0i, [], [], tolerance );
        else
            % Have we reduced the model, then we can rootfind faster
            if ( isfield( ar.config, 'C_rootfinding' ) && ( ar.config.C_rootfinding == 1 ) )
                % Have we compiled with rootfinding capabilities, then use the direct method inside the C-file
                ar.model(jm).(condis)(jc).x0_override = x0i;
                feval(ar.fkt, ar, true, ar.config.useSensis, true, false, condis, threads, 2);
                xnew = ar.model(jm).(condis)(jc).xFineSimu(end,:);
                maxDiff = max(abs(fn(xnew)));
            else
                % Otherwise use the MATLAB implementation
                [xnew, maxDiff] = NewtonRaphson( fn, x0i, [], [], tolerance );
            end
        end
    else
        % Estimate initials in steady state
        if ( debug )
            opts = optimset('TolX', tolerance, 'TolFun', tolerance, 'Jacobian', 'On', 'Display', 'Iter', 'MaxIter', 1e5 ); %, 'DerivativeCheck', 'On'
        else
            opts = optimset('TolX', tolerance, 'TolFun', tolerance, 'Jacobian', 'On', 'Display', 'Off', 'MaxIter', 1e5 );
        end        
        
        xnew = lsqnonlin( fn, x0, x0i, [], opts );
        maxDiff = max(abs(fn(xnew)));
    end
    
    if ( useConserved )
        xnew = totals + map*xnew.';
    end
    xnew = xnew.';
    
    % Calculate sensitivities via implicit function theorem
    if ( useConserved )
        N       = ar.model(jm).N;
        dvdx    = (ar.model(jm).pools.mapping.'*ar.model(jm).(condis)(jc).dvdxNum.').';
        dvdp    = ar.model(jm).(condis)(jc).dvdpNum;
        dfdx    = N(ar.model(jm).pools.usedStates, :) * dvdx; % This is incorrect!
        dfdp    = N(ar.model(jm).pools.usedStates, :) * dvdp;
        S       = - pinv(dfdx) * dfdp;
        S       = ar.model(jm).pools.mapping*S;
    else
        dfdx    = ar.model(jm).(condis)(jc).dfdxNum + 0;
        dfdp    = ar.model(jm).(condis)(jc).dfdpNum + 0;
        S       = pinv(-dfdx)*dfdp;
    end
    
    if ( maxDiff > tolerance )
        warning( 'Failure to converge when rootfinding for model %d, condition %d', jm, jc );
    end
    
    % Remove the override after determination
    ar.model(jm).(condis)(jc).x0_override = [];
    
    if ( debug )        
        totalError = sum(sum((Sref-S).^2)) + sum((xnew-xref).^2);
        disp( 'x found by rootfinding' );
        xnew %#ok
        disp( 'x found by simulating a long time' );
        xref %#ok

        disp( 'S found by rootfinding' );
        S %#ok
        disp( 'Sref found by simulating a long time' );
        Sref %#ok                          
        if ( totalError > ar.config.eq_tol )
            failedCheck = 1;
            
            error( 'Steady state equilibration is different from rootfinding' );      
        else
            disp( 'Test passed' );
        end
    end
end

% dxdts are squared to generate minimum for small dxdt
function [res, J] = merit(x0, jm, jc, condis, threads)
    global ar;
    ar.model(jm).(condis)(jc).x0_override = x0;
    
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, condis, threads, 1);
    res = ar.model(jm).(condis)(jc).dxdt + 0;
    
    if ( nargout > 1 )
        J = ar.model(jm).(condis)(jc).dfdxNum + 0;
    end
end

% dxdts are squared to generate minimum for small dxdt in the presence of
% conservation relations
function res = meritConserved(x0c, jm, jc, map, totals, condis, threads)
    global ar;
    ar.model(jm).(condis)(jc).x0_override = totals + map*x0c.';
    feval(ar.fkt, ar, true, ar.config.useSensis, true, false, condis, threads, 1);
    
    res = ar.model(jm).(condis)(jc).dxdt + 0;
end

function [x, maxDiff] = NewtonRaphson( func, initial, lb, ub, xtol )   
    lastX = inf(size(initial));
    x = initial;
    i = 1;
    while( max( abs( lastX - x ) ) > xtol )
        lastX = x;
        [f,J] = func( x );
        %fprintf( 'Iteration %d; function value: %g\n', i, sum(f.^2) );
        x = x - (pinv(J) * f.').';
    end
    maxDiff = max(abs(f));
end