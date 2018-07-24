% fit sequence
%
% arFits(ps, append)
%
% ps:                           parameter values      
% append:                       [false]
% dynamic_only                  [true]

function arFits(ps, append, dynamic_only, tillConv)
if(~exist('tillConv','var') | isempty(tillConv))
    tillConv = false;
end

global ar

if(~exist('append','var'))
    append = false;
end
if(~exist('dynamic_only','var'))
    dynamic_only = false;
end

n = size(ps,1);

pReset = ar.p;
arCalcMerit
chi2Reset = arGetMerit;

if(append && isfield(ar,'ps') && isfield(ar, 'chi2s') && isfield(ar, 'exitflag') && isfield(ar, 'timing') && isfield(ar, 'fun_evals'))
    ar.ps = [nan(n, length(ar.p)); ar.ps];
    ar.chi2s = [nan(1,n) ar.chi2s];
    ar.timing = [nan(1,n) ar.timing];
    ar.exitflag = [-ones(1,n) ar.exitflag];
    ar.fun_evals = [nan(1,n) ar.fun_evals];
else
    ar.ps = nan(n, length(ar.p));
    ar.chi2s = nan(1,n);
    ar.timing = nan(1,n);
    ar.fun_evals = nan(1,n);
    ar.exitflag = -ones(1,n);
end

ar.ps_start = ps;

if(~dynamic_only)
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
    if(ar.config.fiterrors == 1 || sum(ar.qError==1 & ar.qFit==1)>0)
        text = sprintf('best minimum: -2*log(L) = %f (old = %f)', ...
            min(2*ar.ndata*log(sqrt(2*pi)) + ar.chi2s), 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    else
        text = sprintf('best minimum: chi^2 = %f (old = %f)', min(ar.chi2s), chi2Reset);
    end
    arWaitbar(j, n, text);
    ar.p = ps(j,:);
    tic;
    try
        if(dynamic_only)
            arFitObs(true,2);
        end
        if(tillConv)
            arFitTillConv(false);
        else
%             arFit(true);
            arFitTransient;
        end
        ar.ps(j,:) = ar.p;
        ar.chi2s(j) = ar.chi2fit;
        ar.exitflag(j) = ar.fit.exitflag;
        ar.fun_evals(j) = ar.fevals;
    catch exception
        fprintf('fit #%i: %s\n', j, exception.message);
    end
    ar.timing(j) = toc;
end
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

if(sum(q_select)<6)
    figure(2)
    plotmatrix(ar.ps(:,q_select), 'x');
end

%% Plot
% 
% figure(1)
% 
% [chi2s, isorted] = sort(ar.chi2s);
% exitflag = ar.exitflag(isorted);
% 
% ar.chi2s_sorted = chi2s;
% ar.ps_sorted = ar.ps(isorted,:);
% 
% if(ar.config.fiterrors == 1)    
% %     plot(sort(2*ar.ndata*log(sqrt(2*pi)) + ar.chi2s), 'o--');
%     semilogy(chi2s - ar.chi2fit + 1, '--');
%     hold on
%     h = semilogy(find(exitflag>0), chi2s(exitflag>0) - ar.chi2fit + 1, 'o');
%     plot(xlim, [1 1], 'k--');
%     hold off
%     ylabel('-2*log(L) + const');
% else
%     semilogy(chi2s, '--');
%     hold on
%     h = semilogy(find(exitflag>0), chi2s(exitflag>0), 'o');
%     plot(xlim, [ar.chi2fit ar.chi2fit], 'k--');
%     hold off
%     ylabel('\chi^2');
% end
% xlabel('fits sorted');
% title(sprintf('%i fits in total, %i without errors, %i converged', ...
%     length(exitflag), sum(~isnan(chi2s)) ,sum(exitflag>0)));
% legend(h, 'converged fits');
% 
% % T = clusterdata(log10(chi2s(~isnan(chi2s))'), 1)
% % Y = pdist(log10(chi2s(~isnan(chi2s))'));
% % Z = linkage(Y);
% % figure(3)
% % [H,T] = dendrogram(Z,'colorthreshold','default');
% 
% figure(2)
% subplot(1,2,1)
% hist(ar.timing, 50);
% xlabel('fit time / sec.');
% title(sprintf('total time for %i fits %s', ...
%     length(ar.chi2s), secToHMS(sum(ar.timing(~isnan(ar.timing))))));
% subplot(1,2,2)
% 
% hist(ar.fun_evals, 50);
% xlabel('number of function evaluations');