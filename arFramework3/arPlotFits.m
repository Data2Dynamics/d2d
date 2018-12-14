% arPlotFits(q_select, i_fits)
% 
% Plot results of multistart optimization. This function is intended to be used after arFitLHS or arFits.
% 
% q_select  [ar.qFit == 1]     which parameters are plotted against each other.
% i_fits    [1:size(ar.ps,1)]  which fits are analyzed 
% 
% This function plots three figures visualizing results of multistart
% optimization:
% figure 1: overview about likelihood values before and after the
% multistart optimization. Waterfall plot.
% figure 2: histogram of run time, histogram of function evaluations and
% plot of gradient length.
% figure 3: scatter plots of parameter values (q_select) in the i_fits.

function arPlotFits(q_select, i_fits)

global ar

if(~exist('q_select','var'))
    q_select = ar.qFit==1;
end
if(~exist('i_fits','var'))
    i_fits = 1:size(ar.ps,1);
end

%% chi^2 plots

chi2constr = ar.chi2s + ar.chi2sconstr;
chi2constr_start = ar.chi2s_start + ar.chi2sconstr_start;
[chi2constr, isorted] = sort(chi2constr);
chi2constr_start = chi2constr_start(isorted);

ar.chi2s_sorted = ar.chi2s(isorted);
ar.chi2s_start_sorted = ar.chi2s_start(isorted);
ar.chi2sconstr_sorted = ar.chi2sconstr(isorted);
ar.chi2sconstr_start_sorted = ar.chi2sconstr_start(isorted);

ar.ps_sorted = ar.ps(isorted,:);
exitflag = ar.exitflag(isorted);

chi2min = min([ar.chi2s ar.chi2fit]);
constrmin = min([ar.chi2sconstr ar.chi2constr]);
chi2constrmin = min([ar.chi2s+ar.chi2sconstr ar.chi2fit+ar.chi2constr]);

ar.ps_start_sorted = ar.ps_start(isorted,:);
optim_krit = ar.optim_crit(isorted);

dchi2 = chi2inv(0.95, 1);

figure(1); clf;
pos1 = get(gcf,'Position');
nsub = sum([ar.ndata>0 ar.nconstr>0]);
if(nsub==2)
    nsub = 3;
end
isub = 1;

if(ar.ndata>0)
    subplot(1,nsub,isub);
    
    semilogy(ar.chi2s_sorted - chi2min + 1, '--');
    hold on
    h = semilogy(find(exitflag>0), ar.chi2s_sorted(exitflag>0) - chi2min + 1, 'o');
    h2 = semilogy(ar.chi2s_start_sorted - chi2min + 1, 'x--','Color', [.6 .6 .6]);
    xlim([0 length(chi2constr)+1])
    plot(xlim, [1 1], 'k--');
    plot(xlim, [dchi2 dchi2], 'k:');
    hold off
    title('likelihood');
    xlabel('run index (sorted by likelihood)');
    
    isub = isub + 1;
end

if(ar.nconstr>0)
    subplot(1,nsub,isub);
    
    semilogy(ar.chi2sconstr_sorted - constrmin + 1, '--');
    hold on
    h = semilogy(find(exitflag>0), ar.chi2sconstr_sorted(exitflag>0) - constrmin + 1, 'o');
    h2 = semilogy(ar.chi2sconstr_start_sorted - constrmin + 1, 'x--','Color', [.6 .6 .6]);
    xlim([0 length(chi2constr)+1])
    plot(xlim, [1 1], 'k--');
    plot(xlim, [dchi2 dchi2], 'k:');
    hold off
    title('constraints');
    xlabel('run index (sorted by likelihood)');
    
    isub = isub + 1;
end

if(ar.ndata>0 && ar.nconstr>0)
    subplot(1,nsub,isub);
    
    semilogy(chi2constr - chi2constrmin + 1, '--');
    hold on
    h = semilogy(find(exitflag>0), chi2constr(exitflag>0) - chi2constrmin + 1, 'o');
    h2 = semilogy(chi2constr_start - chi2constrmin + 1, 'x--','Color', [.6 .6 .6]);
    xlim([0 length(chi2constr)+1])
    plot(xlim, [1 1], 'k--');
    plot(xlim, [dchi2 dchi2], 'k:');
    hold off
    title('likelihood + constraints');
    xlabel('run index (sorted by likelihood)');
end
legend([h h2], 'converged fits', 'initial objective function value', 'Location','Best');


figure(2); clf;
subplot(2,2,1);
hist(ar.timing, 50);
xlabel('time / sec / fit');
title({sprintf('run time for %i fits', ...
    sum(~isnan(ar.timing)))
    sprintf('(cumulative time: %s)', ...
    secToHMS(sum(ar.timing(~isnan(ar.timing)))))});

subplot(2,2,2);
hist(ar.fun_evals, 50);
xlabel('number of function evaluations / fit');
title({sprintf('number of function evaluations for %i fits', ...
    sum(~isnan(ar.fun_evals)))
    sprintf('(cumulative number: %i)', ...
    sum(ar.fun_evals(~isnan(ar.fun_evals))))});

subplot(2,2,[3 4]);
semilogy(optim_krit, 'o--');
xlabel('run index (sorted by likelihood)');
ylabel('gradient length');
title('first order optimality criterion');
set(gcf,'Position',get(gcf,'Position')+[pos1(3),0,0,0]); %shift horizontal position by the horizontal width of the 1st figure

%% scatter plots

nji = sum(q_select);
if(nji < 6)
    figure(3);
    
    ji = find(q_select);
    scount = 1;
    for j1 = 1:nji
        for j2 = 1:nji
            if(j1>j2)
                subplot(nji,nji,scount);
                plot(ar.ps_start_sorted(i_fits,ji(j2)), ar.ps_start_sorted(i_fits,ji(j1)), 'ro');
                hold on
                for jn = i_fits
                    plot([ar.ps_start_sorted(jn,ji(j2)) ar.ps_sorted(jn,ji(j2))], ...
                        [ar.ps_start_sorted(jn,ji(j1)) ar.ps_sorted(jn,ji(j1))], '--', 'Color', [.8 .8 .8]);
                    
                end
                plot(ar.ps_sorted(i_fits,ji(j2)), ar.ps_sorted(i_fits,ji(j1)), 'bx');
                hold off
                xlim([ar.lb(ji(j2)) ar.ub(ji(j2))]);
                ylim([ar.lb(ji(j1)) ar.ub(ji(j1))]);
            elseif(j1==j2)
                subplot(nji,nji,scount);
                hist_norm(ar.ps_sorted(i_fits,ji(j1)), linspace(ar.lb(ji(j1)), ar.ub(ji(j1)), 20));
                xlim([ar.lb(ji(j1)) ar.ub(ji(j1))]);
            end
            
            if(scount>nji*(nji-1))
                xlabel(strrep(ar.pLabel{ji(j2)}, '_', '\_'));
            end
            if(j2==1)
                ylabel(strrep(ar.pLabel{ji(j1)}, '_', '\_'));
            end
            scount = scount + 1;
        end
    end
end
