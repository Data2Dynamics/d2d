% plot scan of likelihood

function arPlotScan(jks, savetofile)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(~cellfun(@isempty, ar.scan.ps));
end
if(~exist('savetofile','var'))
    savetofile = false;
end

labelfontsize = 10;
overplot = 0.1;

sumples = 0;
for j=jks
    if(~isempty(ar.scan.ps{j}))
        sumples = sumples + 1;
    end
end

if(~exist('ncols', 'var') || isempty(ncols))
    ncols = ceil(sumples^(0.4))+1;
    nrows = ceil(sumples/ncols);
end
if(~exist('nrows', 'var') || isempty(nrows))
    nrows = ceil(sumples/ncols);
end

h = myRaiseFigure('likelihood scan');
clf
set(h, 'Color', [1 1 1]);

count = 1;

minchi2 = Inf;
for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        minchi2 = min([minchi2 min(ar.scan.chi2s{jk})]);
    end
end
chi2curr = ar.chi2fit;
if(ar.config.fiterrors == 1)
    minchi2 = 2*ar.ndata*log(sqrt(2*pi)) + minchi2;
    chi2curr = 2*ar.ndata*log(sqrt(2*pi)) + chi2curr;
end

dchi2 = chi2inv(1-ar.ppl.alpha_level, ar.ppl.ndof);

for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        g = subplot(nrows,ncols,count);
        
        ps = ar.scan.ps{jk};
        chi2s = ar.scan.chi2s{jk};
        
        if(ar.config.fiterrors == 1)
            chi2s = 2*ar.ndata*log(sqrt(2*pi)) + chi2s;
        end
        
        plot(ps, chi2s, 'k-', 'LineWidth', 1)
        
        ax1 = g;
%         ax2 = axes('Position',get(ax1,'Position'),...
%             'XAxisLocation','top',...
%             'YAxisLocation','right',...
%             'Color','none',...
%             'XColor','r','YColor','r');
        
%         chi2factor = 1.3;
%         ylimmax = minchi2+dchi2*chi2factor;
%         ylim(ax2, [minchi2-(dchi2*chi2factor*overplot) ylimmax+(dchi2*chi2factor*overplot)]);

%         line(ps, chi2s, 'LineWidth', 1, 'Parent', ax2, 'Color', [1 0 0])
        
        % optimum
        line(ar.p(jk), chi2curr, 'Marker', '*', 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerSize', 8, ...
            'Parent', ax1)
        
%         % thresholds
%         line(xlim, [0 0]+minchi2+dchi2, 'LineWidth', 1, 'Parent', ax1, 'Color', [1 0 0], 'LineStyle', '--')
%         if(count == 1)
%             text(mean(xlim), minchi2+dchi2, sprintf('%2i%%', (1-ar.ppl.alpha_level)*100), 'Color', 'r', ...
%                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', labelfontsize, ...
%                 'Parent', ax1)
%         end

        xlabel(ax1, ['log_{10}(' myNameTrafo(ar.pLabel{jk}) ')'])
        spacedAxisLimits(g, overplot)
%         box(g, 'off')
        
        if(mod(count-1,ncols)==0)
            if(ar.config.fiterrors == 1)
                ylabel(ax1, '-2*log(L)');
            else
                ylabel(ax1, '\chi^2');
            end
        end
        
        count = count + 1;
    end
end

% save
if(savetofile && exist(pleGlobals.savePath, 'dir'))
    pleGlobals.figPathMulti{jk} = [pleGlobals.savePath '/multi_plot'];
    save([pleGlobals.savePath '/results.mat'], 'pleGlobals');
    saveas(gcf, [pleGlobals.savePath '/multi_plot'], 'fig')
    print('-depsc2', [pleGlobals.savePath '/multi_plot']);
end


function h = myRaiseFigure(figname)
global pleGlobals
openfigs = get(0,'Children');

figcolor = [1 1 1];

if(isfield(pleGlobals, 'fighandel_multi') && ~isempty(pleGlobals.fighandel_multi) && ...
        pleGlobals.fighandel_multi ~= 0 && ...
        sum(pleGlobals.fighandel_multi==openfigs)>0 && ...
        strcmp(get(pleGlobals.fighandel_multi, 'Name'), figname))
    
    h = pleGlobals.fighandel_multi;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.1 0.6 0.8]);
    set(h,'Color', figcolor);
    pleGlobals.fighandel_multi = h;
end

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');

function spacedAxisLimits(g, overplot)
[xmin xmax ymin ymax] = axisLimits(g);
xrange = xmax - xmin;
if(xrange == 0)
    xrange = 1;
end
yrange = ymax - ymin;
if(yrange == 0)
    yrange = 1;
end
xlim(g, [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
ylim(g, [ymin-(yrange*overplot) ymax+(yrange*overplot)]);

function [xmin xmax ymin ymax] = axisLimits(g)
p = get(g,'Children');
xmin = nan;
xmax = nan;
ymin = nan;
ymax = nan;
for j = 1:length(p)
    if(~strcmp(get(p(j), 'Type'), 'text'))
        xmin = min([xmin toRowVector(get(p(j), 'XData'))]);
        xmax = max([xmax toRowVector(get(p(j), 'XData'))]);
        %         get(p(j), 'UData')
        %         get(p(j), 'LData')
        %         set(p(j), 'LData', get(p(j), 'LData')*2)
        if(strcmp(get(p(j), 'Type'),'hggroup'))
            ymin = min([ymin toRowVector(get(p(j), 'YData'))-toRowVector(get(p(j), 'LData'))]);
            ymax = max([ymax toRowVector(get(p(j), 'YData'))+toRowVector(get(p(j), 'UData'))]);
        else
            ymin = min([ymin toRowVector(get(p(j), 'YData'))]);
            ymax = max([ymax toRowVector(get(p(j), 'YData'))]);
        end
    end
end



function b = toRowVector(a)
b = a(:)';


