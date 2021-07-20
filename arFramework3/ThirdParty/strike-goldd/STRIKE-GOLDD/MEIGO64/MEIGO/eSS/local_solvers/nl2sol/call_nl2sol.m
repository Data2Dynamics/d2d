% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/local_solvers/nl2sol/call_nl2sol.m 1983 2015-02-20 10:55:44Z attila $
%
function [x,fval,exitflag] = call_nl2sol(x0,grad, ndata, lb, ub,opts)

ydata = zeros(ndata,1);

% disp('ndata in call_nl2sol:')
% disp(ndata)

%nl2sol_tuned is a tuned version of that comes with OPTI Toolbox.
% the tuning is essential if gradient is not provided: the perturbation parameter used to estimate the gradient and the Hessian should be adjusted to the precision of the computation of the objective (e.g. the relative tolerance of CVODES.)
% the perturbation parameter

		
[x, fval, exitflag, iter, feval] = nl2sol_v2_1(@objf_nl2sol,grad,x0,ydata,lb,ub,opts); % version by Attila - needs testing
%[x, fval, exitflag, iter, feval] = nl2sol(@objf_nl2sol,grad,x0,ydata,lb,ub,opts);

% be sure:
R = objf_nl2sol(x);
fval = R(:)'*R(:);

switch(exitflag)
    case {3,4}
        info.Status = 'Optimal';
    case 0
        info.Status = 'Exceeded Iterations';
    case -1
        info.Status = 'Infeasible / Could not Converge';
    case -2
        info.Status = 'NL2SOL Error';
    case -5
        info.Status = 'User Exited';
    case 7
        info.Status = 'Singular convergence';
    case 8
        info.Status = 'False convergence';
    otherwise        
        info.Status = 'NL2SOL exited';
        if fval == 0
            fval = inf;
        end
end
fprintf('NL2SOL ended after %d iteration and %d function evaluation:\n', iter,feval);
fprintf('       Results: %s\n',info.Status);
