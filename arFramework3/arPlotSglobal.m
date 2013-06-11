% Plot global sensitivities of residuals
%
% arPlotSRESglobal

function arPlotSglobal(ip, ires)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('ip','var') || isempty(ip))
    ip = 1:size(ar.sres,2);
end

if(~exist('ires','var') || isempty(ires))
    ires = 1:size(ar.sres,1);
end

arChi2(true);

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;
overplot = 0.1;

figure(1)

% rows and cols
np = length(ip);
[nrows, ncols] = NtoColsAndRows(np, rowstocols);

ccount = 1;
for jp = ip
    g = subplot(nrows,ncols,ccount);
    mySubplotStyle(g, labelfontsize, labelfonttype);
    
    plot(g, ar.sres(ires,jp), 'b-o');
    if(isfield(ar,'sresFD'))
        hold(g, 'on');
        plot(g, ar.sresFD(ires,jp), 'r--*');
        hold(g, 'off');
    end
    
    spacedAxisLimits(g, overplot);
    title(g, myNameTrafo(ar.pLabel{jp}));
    if(ccount == 1 && isfield(ar,'sresFD'))
        legend(g, {'SE','FD'});
    end
    
    if(ccount == (nrows-1)*ncols + 1)
        xlabel(g, 'residual');
        ylabel(g, 'sensitivity');
    end
    
    ccount = ccount + 1;
end

if(isfield(ar,'sresFD'))
    figure(2)
    
    mySubplotStyle(gca, labelfontsize, labelfonttype);
    
    semilogy(min(abs(ar.sres(ires,ip) - ar.sresFD(ires,ip))./abs(ar.sres(ires,ip)), ...
        abs(ar.sres(ires,ip) - ar.sresFD(ires,ip))), 'x-')
    legend(myNameTrafo(ar.pLabel(ip)))
end

% % constrains
% if(isfield(ar,'constr'))
%     figure(3)
%     
%     % rows and cols
%     np = size(ar.sconstr,2);
%     [nrows, ncols] = NtoColsAndRows(np, rowstocols);
%     
%     
%     for jp = 1:np
%         g = subplot(nrows,ncols,jp);
%         mySubplotStyle(g, labelfontsize, labelfonttype);
%         
%         plot(g, ar.sconstr(:,jp), 'b-o');
%         if(isfield(ar,'sconstrFD'))
%             hold(g, 'on');
%             plot(g, ar.sconstrFD(:,jp), 'r--*');
%             hold(g, 'off');
%         end
%         
%         spacedAxisLimits(g, overplot);
%         title(g, myNameTrafo(ar.pLabel{jp}));
%         if(jp == 1 && isfield(ar,'sconstrFD'))
%             legend(g, {'SE','FD'});
%         end
%         
%         if(jp == (nrows-1)*ncols + 1)
%             xlabel(g, 'constraint');
%             ylabel(g, 'sensitivity');
%         end
%     end
%     
%     if(isfield(ar,'sconstrFD'))
%         figure(4)
%         
%         mySubplotStyle(gca, labelfontsize, labelfonttype);
%         
%         semilogy(min(abs(ar.sconstr - ar.sconstrFD)./abs(ar.sconstr), ...
%             abs(ar.sconstr - ar.sconstrFD)), 'x-')
%         legend(myNameTrafo(ar.pLabel))
%     end
% end

%% sub-functions

function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);


function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);



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


