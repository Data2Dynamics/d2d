% modified fit sequence
%
% arFits_modified(ps, append)
%
% ps:                           parameter values      
% append:                       [false]
% dynamic_only                  [true]

function arFits_modified(ps, append, dynamic_only, plot_summary)

global ar

if(~exist('append','var'))
    append = false;
end
if(~exist('dynamic_only','var'))
    dynamic_only = true;
end
if(~exist('plot_summary','var'))
    plot_summary = false;
end

ar.config.useFitErrorCorrection = false;

n = size(ps,1);

arChi2(true);
pReset = ar.p;
chi2Reset = ar.chi2fit;

if(append && isfield(ar,'ps') && isfield(ar, 'chi2s') && ...
        isfield(ar, 'exitflag') && isfield(ar, 'timing') && ...
        isfield(ar, 'fun_evals'))
    ar.ps = [nan(n, length(ar.p)); ar.ps];
    ar.ps_start = [ps; ar.ps_start];
    ar.chi2s = [nan(1,n) ar.chi2s];
    ar.chi2s_start = [nan(1,n) ar.chi2s_start];
    ar.timing = [nan(1,n) ar.timing];
    ar.exitflag = [-ones(1,n) ar.exitflag];
    ar.fun_evals = [nan(1,n) ar.fun_evals];
    ar.optim_crit = [nan(1,n) ar.optim_crit];
else
    ar.ps = nan(n, length(ar.p));
    ar.ps_start = ps;
    ar.chi2s = nan(1,n);
    ar.chi2s_start = nan(1,n);
    ar.timing = nan(1,n);
    ar.fun_evals = nan(1,n);
    ar.optim_crit = nan(1,n);
    ar.exitflag = -ones(1,n);
end

if(dynamic_only)
    q_select = ar.qFit==1 & ar.qDynamic==1;
else
    q_select = ar.qFit==1;
end

if(sum(q_select)<6)
    figure(1)
    plotmatrix(ps(:,q_select), 'x');
end

arWaitbar(0);
for j=1:n
    if(ar.config.fiterrors == 1)  
        text = sprintf('best minimum: -2*log(L) = %f (old = %f)', ...
            min(2*ar.ndata*log(sqrt(2*pi)) + ar.chi2s), 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    else
        text = sprintf('best minimum: chi^2 = %f (old = %f)', min(ar.chi2s), chi2Reset);
    end
    arWaitbar(j, n, text);
    ar.p = ps(j,:);
    tic;

    try
        arChi2(true);
        ar.chi2s_start(j) = ar.chi2fit;
        i_max = 2;
        ar.qFit(:) = 1;

	%Adjust borders, error-parameters / non-error-parameters
        for i = 1:length(ar.p)
            if(ar.qError(i) == 0)
                ar.ub(i) = +3;
                ar.lb(i) = -5;
            elseif(ar.qError(i) == 1)
                ar.ub(i) = +2;
                ar.lb(i) = -2;
            end
        end

        for i = 1:i_max

            %Choose optimizer
            %ar.config.optimizer = 1;

            ar.config.rtol = 1e-6;
            ar.config.atol = 1e-6;
            ar.config.optim.TolX = 0;
            ar.config.optim.TolFun = 0;

            if(i == 1 || i == 2)
                ar.config.optim.MaxIter = 10;
                for nn = 1:100
                    if(i==1 && nn==1)
                        for l = 1:length(ar.p)
                            if(ar.qDynamic(l) == 0)
                                ar.qFit(l) = 0;
                            end
                        end
			fprintf('Run %i, non-dynamic parameters = 0\n',i);
                    elseif(i==2 && nn==1)
                        for l = 1:length(ar.p)
                            if(ar.qDynamic(l) == 0)
                                ar.qFit(l) = 1;
                            end
                        end
			fprintf('Run %i, non-dynamic parameters = 1\n',i);
                    end;

                    if(rand >= 0.75)
                        fprintf('All Parameters are released\n');
                        for k = 1:length(ar.p)
                            if(ar.qDynamic(k) == 1)
                                ar.qFit(k) = 1;
                            end;
                        end;
                        if(i==1)
                            for l = 1:length(ar.p)
                                if(ar.qDynamic(l) == 0)
                                    ar.qFit(l) = 0;
                                end
                            end
                        elseif(i==2)
                            for l = 1:length(ar.p)
                                if(ar.qDynamic(l) == 0)
                                    ar.qFit(l) = 1;
                                end
                            end
                        end;
                    end;

                    tmp = ar.chi2fit;

                    try
                        arFit();
                    catch exception2
                        %Try to save crashed run by disturbing last parameter-set
                        fprintf('Breakdown catched\n');
                        ar.qFit(:) = 1;
                        arDisturb();
                        arFit();
                    end

                    ar.p_tmp(n,:) = ar.p;

                    if(ar.chi2fit-tmp == 0)
                        break;
                    end

                    for k = 1:length(ar.p)
                        if((ar.p(k) < ar.lb(k)*0.9 || ar.p(k) > ar.ub(k)*0.9) && ar.qFit(k) == 1)
                            if(ar.p(k)>=1)
                                ar.qFit(k) = 0;
                                fprintf('Para %i is fixed(high) \n ',k);
                            elseif(ar.p(k)<=-1)
                                ar.qFit(k) = 0;
                                fprintf('Para %i is fixed(low) \n ',k);
                            end;    
                        end;
                    end;
                end
            end
        end

        fprintf('Final fit..\n');
        ar.config.optim.MaxIter = 10;
        ar.qFit(:) = 1;
        for p = 1:5
            tmp = ar.chi2fit;
            arFit();
             if(ar.chi2fit-tmp == 0)
                 break;
             end
        end
        fprintf('\nFinished..\n');

        ar.ps(j,:) = ar.p;
        ar.chi2s(j) = ar.chi2fit;
        ar.exitflag(j) = ar.fit.exitflag;
        ar.fun_evals(j) = ar.fevals;
        ar.optim_crit(j) = ar.fit.output.firstorderopt;
    catch exception
        fprintf('fit #%i: %s\n', j, exception.message);
    end
    ar.timing(j) = toc;
end
fprintf('total fitting time: %fsec\n', sum(ar.timing(~isnan(ar.timing))));
fprintf('mean fitting time: %fsec\n', 10^mean(log10(ar.timing(~isnan(ar.timing)))));
arWaitbar(-1);

if(chi2Reset>min(ar.chi2s))
    [chi2min,imin] = min(ar.chi2s);
    ar.p = ar.ps(imin,:);
    if(ar.config.fiterrors == 1)
        fprintf('selected best fit #%i with -2*log(L) = %f (old = %f)\n', ...
            imin, 2*ar.ndata*log(sqrt(2*pi)) + chi2min, 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    else
        fprintf('selected best fit #%i with chi^2 = %f (old = %f)\n', imin, chi2min, chi2Reset);
    end
else
    fprintf('did not find better fit\n');
    ar.p = pReset;
end
arChi2(false);

if plot_summary
    if(sum(q_select)<6)
        figure(2)
        plotmatrix(ar.ps(:,q_select), 'x');
    end
    
    arPlotFits;
end
