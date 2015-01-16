% Plot models sensitivities of residuals
%
% arPlotSY

function arPlotSRES

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

% arChi2(true);

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;
overplot = 0.1;

fcount = 1;
for jm = 1:length(ar.model)
    nd = length(ar.model(jm).data);
    for jd = 1:nd
        myRaiseFigure(jm, ['SRES: ' ar.model(jm).data(jd).name ' - ' ar.model(jm).data(jd).checkstr], fcount, jd);
        
        % rows and cols
        [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols);
        
        np = length(ar.model(jm).data(jd).p);
        for jy = 1:ny
            g = subplot(nrows,ncols,jy);
            arSubplotStyle(g, labelfontsize, labelfonttype);
            
            legendhandle = zeros(1,np);
            
            for jp = 1:np
                linestyle = myLineStyle(np,jp);
                ltmp = plot(g, ar.model(jm).data(jd).tExp, ar.model(jm).data(jd).sres(:,jy,jp), linestyle{:}, 'Marker', 'o');
                legendhandle(jp) = ltmp;
                hold(g, 'on');
                if(isfield(ar.model(jm).data(jd), 'sresFD'))
                    plot(g, ar.model(jm).data(jd).tExp, ar.model(jm).data(jd).sresFD(:,jy,jp), linestyle{:}, 'Marker', '*');
                end
            end
            hold(g, 'off');
            
            spacedAxisLimits(g, overplot);
            title(g, myNameTrafo(ar.model(jm).data(jd).y{jy}));
            if(jy == 1)
                legend(g, legendhandle, myNameTrafo(ar.model(jm).data(jd).p));
            end
            
            if(jy == (nrows-1)*ncols + 1)
                xlabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).tUnits{3}, ar.model(jm).data(jd).tUnits{2}));
                ylabel(g, 'sensitivity');
            end
        end      
        fcount = fcount + 1;
    end
end

if(ar.config.fiterrors == 1)
    arPlotSRESERR
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




function h = myRaiseFigure(m, figname, jk, jd)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

if(isfield(ar.model(m).data(jd), 'fighandel_sres') && ~isempty(ar.model(m).data(jd).fighandel_sres) && ...
        ar.model(m).data(jd).fighandel_sres ~= 0 && sum(ar.model(m).data(jd).fighandel_sres==openfigs)>0 && ...
        strcmp(get(ar.model(m).data(jd).fighandel_sres, 'Name'), figname))
    h = ar.model(m).data(jd).fighandel_sres;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.05+((jk-1)*figdist) 0.45-((jk-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.model(m).data(jd).fighandel_sres = h;
end



function str = myNameTrafo(str)
str = strrep(str, '_', '\_');





function [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols)
global ar
ny = size(ar.model(jm).data(jd).y, 2);
[nrows, ncols] = NtoColsAndRows(ny, rowstocols);



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


