% Plot residuals
%
% arPlotResiduals(saveToFile)
%
% saveToFile    [false]

function arPlotResiduals(saveToFile)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('saveToFile','var'))
    saveToFile = false;
end

h = arRaiseFigure(ar.plot, 'fighandel_res', 'Residuals');

tmpres = [];
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        for y=1:size(ar.model(m).data(d).res, 2)
            tmpres1 = ar.model(m).data(d).res(:,y);
            tmpres = [tmpres; tmpres1(~isnan(tmpres1))]; %#ok<AGROW>
        end
    end
end

g = subplot(1,2,1);
arSubplotStyle(g);
qqplot(tmpres);
box(g, 'on')
axis(g,'equal')
axis(g,'square')
title(g, 'QQ-plot of sample data versus standard normal');
ylabel(g, 'Quantiles of Input Sample');
xlabel(g, 'Standard Normal Quantiles');

g = subplot(1,2,2);
arSubplotStyle(g);
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        for y=1:size(ar.model(m).data(d).res, 2)
            tmpres = ar.model(m).data(d).res(:,y);
            ttmp = ar.model(m).data(d).tExp(~isnan(tmpres));
            tmpres = tmpres(~isnan(tmpres));
			tunique = unique(ttmp); %R2013a compatible
            for jt = 2:length(tunique);
                res1 = tmpres(ttmp==tunique(jt));
                res2 = tmpres(ttmp==tunique(jt-1));
                for jr=1:length(res1)
                    plot(res1(jr), res2, 'b+')
                    hold(g, 'on');
                end
            end
        end
    end
end
plot(xlim*0, ylim, 'r-.');
plot(xlim, ylim*0, 'r-.');
hold(g, 'off');
box(g, 'on')
axis(g,'equal')
axis(g,'square')
title(g, 'Auto-correlation of normalized residuals')
ylabel(g, 'residuals t_i')
xlabel(g, 'residuals t_{i+1}')

if(saveToFile)
    ar.plot.savePath_ResFig = arSaveFigure(h, 'residuals', '/Figures');
end

