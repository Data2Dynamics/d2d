% Plot compact profile Likelihoods
%
% plePlotMultiLLH(ncols, nrows, savetofile)
%
% savetofile        save plot                                       [false]

function plePlotMultiLLH(ncols, nrows, savetofile)

global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end
if(~exist('savetofile','var'))
    savetofile = false;
end

labelfontsize = 10;

sumples = 0;
for j=1:length(pleGlobals.ps)
    if(~isempty(pleGlobals.ps{j}))
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

if(pleGlobals.plot_point && ~pleGlobals.plot_simu)
    strtitle = sprintf('profile likelihood (point-wise)');
elseif(~pleGlobals.plot_point && pleGlobals.plot_simu)
    strtitle = sprintf('profile likelihood (simultaneous)');
else
    strtitle = sprintf('profile likelihood');
end

h = myRaiseFigure(strtitle);
set(h, 'Color', [1 1 1]);

count = 1;
for jk=1:length(pleGlobals.ps)
    if(~isempty(pleGlobals.ps{jk}))
        g = subplot(nrows,ncols,count);
        set(g, 'FontSize', labelfontsize);
        set(g, 'FontName', 'TimesNewRoman');

        % profile
        plot(pleGlobals.ps{jk}(:,jk), transformFromLog(pleGlobals.chi2s{jk}), 'k', 'LineWidth', 1)
        hold on
        
        % thresholds
        if(pleGlobals.plot_point)
            plot(xlim, transformFromLog([0 0]+min(pleGlobals.chi2s{jk})+pleGlobals.dchi2_point), 'r--', 'LineWidth', 1)
        end
        if(pleGlobals.plot_simu)
            plot(xlim, transformFromLog([0 0]+min(pleGlobals.chi2s{jk})+pleGlobals.dchi2), 'r--', 'LineWidth', 1)
        end
        
        % hessian CI
%         sampletmp = linspace(pleGlobals.conf_lb(1,count), pleGlobals.conf_ub(1,count), 100);
%         assym_chi2 = min(pleGlobals.chi2s{count}) + (1/pleGlobals.pstd(count)^2)*(pleGlobals.p(count)-sampletmp).^2;
%         if(isreal(sampletmp))
%             plot(sampletmp, transformFromLog(assym_chi2), '-', 'Color', [.5 .5 .5], 'LineWidth', 1)
%         end
        
        % optimum
        plot(pleGlobals.p(jk), transformFromLog(pleGlobals.chi2), '*', 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerSize', labelfontsize)
        hold off
        
        xlimtmp2 = (max(pleGlobals.ps{jk}(:,jk))-min(pleGlobals.ps{jk}(:,jk)))*0.05;
        if(xlimtmp2>0)
            xlimtmp = [min(pleGlobals.ps{jk}(:,jk))-xlimtmp2 max(pleGlobals.ps{jk}(:,jk))+xlimtmp2];
            xlim(xlimtmp);
        end  
        xlabel(['log_{10}(' myNameTrafo(pleGlobals.p_labels{jk})])
        
        dchi2 = pleGlobals.dchi2_point;
        if(pleGlobals.plot_simu)
            dchi2 = pleGlobals.dchi2;
        end
        
        ylimmax = max(transformFromLog(pleGlobals.chi2s{jk}))*1.1;
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
if(savetofile && exist(pleGlobals.savePath, 'dir'))
    saveas(gcf, [pleGlobals.savePath '/multi_plot_llh'], 'fig')
    print('-depsc2', [pleGlobals.savePath '/multi_plot_llh']);
end


function b = transformFromLog(a)
b = exp(-0.5*a);

function h = myRaiseFigure(figname)
global pleGlobals
openfigs = get(0,'Children');

figcolor = [1 1 1];

if(isfield(pleGlobals, 'fighandel_multi_llh') && ~isempty(pleGlobals.fighandel_multi_llh) && ...
    pleGlobals.fighandel_multi_llh ~= 0 && ...
    sum(pleGlobals.fighandel_multi_llh==openfigs)>0 && ...
    strcmp(get(pleGlobals.fighandel_multi_llh, 'Name'), figname))

    h = pleGlobals.fighandel_multi_llh;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.1 0.6 0.8]);
    set(h,'Color', figcolor);
    pleGlobals.fighandel_multi_llh = h;
end

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');
