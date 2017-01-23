% Plot compact profile Likelihoods
%
% plePlotMultiLLH(ncols, nrows, savetofile)
%
% savetofile        save plot                                       [false]

function plePlotMultiLLH(ncols, nrows, savetofile)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end
if(~exist('savetofile','var'))
    savetofile = false;
end

sumples = 0;
length(ar.ple.ps)
for j=1:length(ar.ple.ps)
    if(~isempty(ar.ple.ps{j}))
        sumples = sumples + 1;
    end
end

if(~exist('ncols', 'var'))
    ncols = ceil(sumples^(0.4))+1;
    nrows = ceil(sumples/ncols);
end
if(~exist('nrows', 'var'))
    nrows = ceil(sumples/ncols);
end

if(ar.ple.plot_point && ~ar.ple.plot_simu)
    strtitle = sprintf('profile likelihood (point-wise)');
elseif(~ar.ple.plot_point && ar.ple.plot_simu)
    strtitle = sprintf('profile likelihood (simultaneous)');
else
    strtitle = sprintf('profile likelihood');
end

h = myRaiseFigure(strtitle);
set(h, 'Color', [1 1 1]);

count = 1;
for jk=1:length(ar.ple.ps)
    if(~isempty(ar.ple.ps{jk}))
        g = subplot(nrows,ncols,count);
        arSubplotStyle(g);

        % profile
        plot(ar.ple.ps{jk}(:,jk), transformFromLog(ar.ple.chi2s{jk}), 'k', 'LineWidth', 1)
        hold on
        
        % thresholds
        if(ar.ple.plot_point)
            plot(xlim, transformFromLog([0 0]+min(ar.ple.chi2s{jk})+ar.ple.dchi2_point), 'r--', 'LineWidth', 1)
        end
        if(ar.ple.plot_simu)
            plot(xlim, transformFromLog([0 0]+min(ar.ple.chi2s{jk})+ar.ple.dchi2), 'r--', 'LineWidth', 1)
        end
        
        % hessian CI
%         sampletmp = linspace(ar.ple.conf_lb(1,count), ar.ple.conf_ub(1,count), 100);
%         assym_chi2 = min(ar.ple.chi2s{count}) + (1/ar.ple.pstd(count)^2)*(ar.ple.p(count)-sampletmp).^2;
%         if(isreal(sampletmp))
%             plot(sampletmp, transformFromLog(assym_chi2), '-', 'Color', [.5 .5 .5], 'LineWidth', 1)
%         end
        
        % optimum
        plot(ar.ple.p(jk), transformFromLog(ar.ple.merit), '*', 'Color', [.5 .5 .5], 'LineWidth', 1)
        hold off
        
        xlimtmp2 = (max(ar.ple.ps{jk}(:,jk))-min(ar.ple.ps{jk}(:,jk)))*0.05;
        if(xlimtmp2>0)
            xlimtmp = [min(ar.ple.ps{jk}(:,jk))-xlimtmp2 max(ar.ple.ps{jk}(:,jk))+xlimtmp2];
            xlim(xlimtmp);
        end  
        xlabel(['log_{10}(' arNameTrafo(ar.ple.p_labels{jk})])
        
        dchi2 = ar.ple.dchi2_point;
        if(ar.ple.plot_simu)
            dchi2 = ar.ple.dchi2;
        end
        
        ylimmax = max(transformFromLog(ar.ple.chi2s{jk}))*1.1;
        ylim([0 ylimmax]);
        
        if(mod(count-1,ncols)==0)
            ylabel('\chi^2_{PL}')
        else
            ylabel('')
            set(gca, 'YTickLabel', {})
        end
        
        count = count + 1;
    end
end

% save
if(savetofile && exist(ar.ple.savePath, 'dir'))
    saveas(gcf, [ar.ple.savePath '/multi_plot_llh'], 'fig')
    print('-depsc2', [ar.ple.savePath '/multi_plot_llh']);
end


function b = transformFromLog(a)
b = exp(-0.5*a);

function h = myRaiseFigure(figname)

openfigs = get(0,'Children');

figcolor = [1 1 1];

if(isfield(ar.ple, 'fighandel_multi_llh') && ~isempty(ar.ple.fighandel_multi_llh) && ...
    ar.ple.fighandel_multi_llh ~= 0 && ...
    sum(ar.ple.fighandel_multi_llh==openfigs)>0 && ...
    strcmp(get(ar.ple.fighandel_multi_llh, 'Name'), figname))

    h = ar.ple.fighandel_multi_llh;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.1 0.6 0.8]);
    set(h,'Color', figcolor);
    ar.ple.fighandel_multi_llh = h;
end
