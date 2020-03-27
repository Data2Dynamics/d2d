% arPlotResiduals(saveToFile)
%
% Plot residuals as quantile-quantile plot and
% autocorrelation (res_{i} versus res_{i+1})
%
%   saveToFile    [false]

function arPlotResiduals(saveToFile)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('saveToFile','var'))
    saveToFile = false;
end

h = myRaiseFigure('Residuals');

tmpres = [];
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        for y=1:size(ar.model(m).data(d).res, 2)
            if ar.model(m).data(d).qFit(y) == 0
                continue
            end
            tmpres1 = ar.model(m).data(d).res(~isnan(ar.model(m).data(d).yExp(:,y)),y);
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
            if ar.model(m).data(d).qFit(y) == 0
                continue
            end
            tmpres = ar.model(m).data(d).res(~isnan(ar.model(m).data(d).yExp(:,y)),y);
            ttmp = ar.model(m).data(d).tExp(~isnan(ar.model(m).data(d).yExp(:,y)));
            ttmp = ttmp(~isnan(tmpres));
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

%% sub-functions



function h = myRaiseFigure(figname)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];

ar.plot.time = now;

if(isfield(ar.plot, 'fighandel_res') && ~isempty(ar.plot.fighandel_res) && ...
        ar.plot.fighandel_res ~= 0 && ...
        sum(ar.plot.fighandel_res==openfigs)>0 && ...
        strcmp(get(ar.plot.fighandel_res, 'Name'), figname))
    
    h = ar.plot.fighandel_res;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.71 0.7 0.25 0.2]);
    set(h,'Color', figcolor);
    ar.plot.fighandel_res = h;
end



