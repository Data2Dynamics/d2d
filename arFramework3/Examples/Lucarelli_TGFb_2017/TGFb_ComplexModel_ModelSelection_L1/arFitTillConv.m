%% Fit till convergence

function arFitTillConv(show_waitbar)
if(~exist('show_waitbar','var') | isempty(show_waitbar))
    show_waitbar = true;
end
global ar
global which_trdog

if(show_waitbar)
    arWaitbar(0);
end
exitflag = 0;
maxIters = 100;

for k = 1:maxIters
    which_trdog = 1;
    try
        arFit(true);
    catch exception
        fprintf(exception.message);
    end
    exitflag = ar.fit.exitflag;

    which_trdog = 0;
    try
        arFit(true);
    catch exception
        fprintf(exception.message);
    end
    exitflag0 = ar.fit.exitflag;

    if exitflag == 2 & exitflag0 ==2
        break
    end
    if(show_waitbar)
        arWaitbar(k,maxIters,sprintf('Current -2*log(L): %f',2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit));
    end
end

if k == maxIters
    disp('Warning: Maximal iterations reached ...');
end
if show_waitbar
    arWaitbar(-1);
end
