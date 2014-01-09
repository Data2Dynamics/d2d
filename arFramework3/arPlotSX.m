% Plot X models sensitivities
%
% arPlotSX

function arPlotSX

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;
overplot = 0.1;

fcount = 1;
for jm = 1:length(ar.model)
    nc = length(ar.model(jm).condition);
    for jc = 1:nc
        myRaiseFigure(jm, ['SX: ' ar.model(jm).name ' - ' ar.model(jm).condition(jc).checkstr], fcount);
        
        % rows and cols
        [ncols, nrows, nu, nx, nz] = myColsAndRows(jm, rowstocols);
        
        np = length(ar.model(jm).condition(jc).p);
        for ju = 1:nu
            g = subplot(nrows,ncols,ju);
            mySubplotStyle(g, labelfontsize, labelfonttype);
            
            legendhandle = zeros(1,np);
            
            for jp = 1:np
                linestyle = myLineStyle(np,jp);
                ltmp = plot(g, ar.model(jm).condition(jc).tFine, ar.model(jm).condition(jc).suFineSimu(:,ju,jp), linestyle{:});
                legendhandle(jp) = ltmp;
                hold(g, 'on');
                % plot(g, ar.model(jm).condition(jc).tExp, squeeze(ar.model(jm).condition(jc).suExpSimu(:,ju,jp)), 'o');
                if(isfield(ar.model(jm).condition(jc), 'suExpSimuFD'))
                    plot(g, ar.model(jm).condition(jc).tExp, ar.model(jm).condition(jc).suExpSimuFD(:,ju,jp), linestyle{2:3}, 'Marker', '*');
                end
            end
            hold(g, 'off');
            
            spacedAxisLimits(g, overplot);
            title(g, myNameTrafo(ar.model(jm).u{ju}));
            if(ju == 1)
                legend(g, legendhandle, myNameTrafo(ar.model(jm).condition(jc).p));
            end
        end
        for jx = 1:nx
            g = subplot(nrows,ncols,jx+nu);
            mySubplotStyle(g, labelfontsize, labelfonttype);
            
            for jp = 1:np
                linestyle = myLineStyle(np,jp);
                plot(g, ar.model(jm).condition(jc).tFine, squeeze(ar.model(jm).condition(jc).sxFineSimu(:,jx,jp)), linestyle{:});
                hold(g, 'on');
                % plot(g, ar.model(jm).condition(jc).tExp, squeeze(ar.model(jm).condition(jc).sxExpSimu(:,jx,jp)), 'o');
                if(isfield(ar.model(jm).condition(jc), 'sxExpSimuFD'))
                    plot(g, ar.model(jm).condition(jc).tExp, ar.model(jm).condition(jc).sxExpSimuFD(:,jx,jp), linestyle{:}, 'Marker', '*');
                end
            end
            hold(g, 'off');
            
            spacedAxisLimits(g, overplot);
            title(g, myNameTrafo(ar.model(jm).x{jx}));
            
            if(jx+nu == (nrows-1)*ncols + 1)
                xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                ylabel(g, 'sensitivity');
            end
        end
        for jz = 1:nz
            g = subplot(nrows,ncols,jz+nu+nx);
            mySubplotStyle(g, labelfontsize, labelfonttype);
            
            for jp = 1:np
                linestyle = myLineStyle(np,jp);
                plot(g, ar.model(jm).condition(jc).tFine, squeeze(ar.model(jm).condition(jc).szFineSimu(:,jz,jp)), linestyle{:});
                hold(g, 'on');
                % plot(g, ar.model(jm).condition(jc).tExp, squeeze(ar.model(jm).condition(jc).szExpSimu(:,jz,jp)), 'o');
                if(isfield(ar.model(jm).condition(jc), 'szExpSimuFD'))
                    plot(g, ar.model(jm).condition(jc).tExp, ar.model(jm).condition(jc).szExpSimuFD(:,jz,jp), linestyle{:}, 'Marker', '*');
                end
            end
            hold(g, 'off');
            
            spacedAxisLimits(g, overplot);
            title(g, myNameTrafo(ar.model(jm).z{jz}));
            
            if(jz+nu+nx == (nrows-1)*ncols + 1)
                xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                ylabel(g, 'sensitivity');
            end
        end
        fcount = fcount + 1;
    end
end



%% sub-functions



function C = myLineStyle(n, j)
farben = lines(n);
zeichen = {':', '-', '--', '-.'};
zeichenindex = mod(floor((j-1)/7)+1, 4)+1;
C = cell(1,3);
C{1} = [zeichen{zeichenindex}];
C{2} = 'Color';
C{3} = farben(j,:);




function h = myRaiseFigure(m, figname, jf)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plots(jf).time = now;

if(isfield(ar.model(m).plots(jf), 'fighandel_sx') && ~isempty(ar.model(m).plots(jf).fighandel_sx) && ...
        ar.model(m).plots(jf).fighandel_sx ~= 0 && sum(ar.model(m).plots(jf).fighandel_sx==openfigs)>0 && ...
        strcmp(get(ar.model(m).plots(jf).fighandel_sx, 'Name'), figname))
    h = ar.model(m).plots(jf).fighandel_sx;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.45+((jf-1)*figdist) 0.45-((jf-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.model(m).plots(jf).fighandel_sx = h;
end



function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);



function [ncols, nrows, nu, nx, nz] = myColsAndRows(jm, rowstocols)
global ar
nu = size(ar.model(jm).u, 2);
nx = size(ar.model(jm).x, 2);
nz = size(ar.model(jm).z, 2);
[nrows, ncols] = NtoColsAndRows(nu+nx+nz, rowstocols);



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);



function spacedAxisLimits(g, overplot)
[xmin, xmax, ymin, ymax] = axisLimits(g);
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



function [xmin, xmax, ymin, ymax] = axisLimits(g)
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


