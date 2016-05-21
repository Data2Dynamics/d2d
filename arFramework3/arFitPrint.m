function arFitPrint

global ar

outputstr = '';
switch ar.fit.exitflag
    case 1
        outputstr = 'Converged to a solution';
    case 2
        outputstr = 'Change in X too small';
    case 3
        outputstr = 'Change in RESNORM too small';
    case 4
        outputstr = 'Computed search direction too small';
    case 0
        outputstr = 'Too many function evaluations or iterations';
    case -1
        outputstr = 'Stopped by output/plot function';
    case -2
        outputstr = 'Bounds are inconsistent';
    case -3
        outputstr = 'Regularization parameter too large (Levenberg-Marquardt)';
    case -4
        outputstr = 'Line search failed';
    case -98
        outputstr = sprintf('Multiple Shooting: constraint strength not controlable > %e\n', 1e20);
    case -99
        outputstr = sprintf('Multiple Shooting: mean constraint violation > %e\n', ar.ms_treshold);
    case 50
        outputstr = 'MaxIter Reached';
    case 51
        outputstr = 'FitCount Reached';
    case 60
        outputstr = 'Trust-Region Radius is Less than Floating-Point Relative Accuracy.';
end

% Special output for Ceres optimizer, since it creates own exitstring (not
% numeric exitflag)
if(ar.config.optimizer==10)
    arFprintf(1, '%s finished after %i iterations:\nExit message: %s\ntotal improvement = %g\n', ...
        ar.config.optimizers{ar.config.optimizer}, ...
        ar.fit.output.iterations, ar.fit.ceresexitmessage, ar.fit.improve);
else    
    arFprintf(1, '%s finished after %i iterations: %s, total improvement = %g\n', ...
        ar.config.optimizers{ar.config.optimizer}, ...
        ar.fit.output.iterations, outputstr, ar.fit.improve);
end
