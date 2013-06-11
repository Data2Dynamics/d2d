% Plot start and end of fits

function arPlotFits(q_select, i_fits)

global ar

if(~exist('q_select','var'))
    q_select = ar.qFit==1;
end
if(~exist('i_fits','var'))
    i_fits = 1:size(ar.ps,1);
end

%% chi^2 plots

[chi2s, isorted] = sort(ar.chi2s);
exitflag = ar.exitflag(isorted);

ar.chi2s_sorted = chi2s;
ar.ps_sorted = ar.ps(isorted,:);

chi2min = min(chi2s);

if(isfield(ar,'chi2s_start'))
    ar.chi2s_start_sorted = ar.chi2s_start(isorted);
end
if(isfield(ar,'ps_start'))
    ar.ps_start_sorted = ar.ps_start(isorted,:);
end

figure(1)

subplot(3,1,[1 2]);
semilogy(chi2s - chi2min + 1, '--');
hold on
h = semilogy(find(exitflag>0), chi2s(exitflag>0) - chi2min + 1, 'o');
if(isfield(ar,'chi2s_start'))
    h2 = semilogy(ar.chi2s_start_sorted - chi2min + 1, 'x--','Color', [.6 .6 .6]);
else
    h2 = [];
end
plot(xlim, [1 1], 'k--');
hold off
if(ar.config.fiterrors == 1)    
    ylabel('-2*log(L) + const');
else
    ylabel('\chi^2');
end
xlabel('fits sorted');
xlim([0 length(chi2s)+1])
title(sprintf('%i fits in total, %i without errors, %i converged', ...
    length(exitflag), sum(~isnan(chi2s)) ,sum(exitflag>0)));
legend([h h2], 'converged fits', 'initial value', 'Location','Best');

subplot(3,2,5);
hist(ar.timing, 50);
xlabel('fit time / sec.');
title(sprintf('total time for %i fits %s', ...
    length(ar.chi2s), secToHMS(sum(ar.timing(~isnan(ar.timing))))));

subplot(3,2,6);
hist(ar.fun_evals, 50);
xlabel('number of function evaluations');

%% scatter plots

nji = sum(q_select);
if(nji < 6)
    figure(2);
    
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
