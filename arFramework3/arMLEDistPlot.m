% check distribution of MLE
%
% arMLEDistPlot(plotChi2dist, saveToFile)
%
% plotChi2dist      [true]
% saveToFile        [false]

function arMLEDistPlot(plotChi2dist, saveToFile)

global ar
global mledist

if(~exist('plotChi2dist','var'))
    plotChi2dist = true;
end
if(~exist('saveToFile','var'))
    saveToFile = false;
end

% constants
labelfontsize = 8;
labelfonttype = 'TimesNewRoman';

h = myRaiseFigure;

timesNdata = 4;

nbins = 50;
subplot(3,3,1)
if(mledist.fiterrors==1)
	xbins = linspace(min(mledist.chi2),max(mledist.chi2), nbins);
else
	xbins = linspace(1,mledist.ndata*timesNdata, nbins);
end
xbinsdiff = diff(xbins);
[nc] = histc(mledist.chi2, xbins);
nc = nc / sum(nc) / xbinsdiff(1);
bar(xbins, nc, 'EdgeColor', 'k', 'FaceColor','w');
arSubplotStyle(gca, labelfontsize, labelfonttype);
if(plotChi2dist)
	if(mledist.fiterrors==1)
		chi2tmp = linspace(min(mledist.chi2),max(mledist.chi2), 100);
	else
		chi2tmp = linspace(0,mledist.ndata*timesNdata, 100);
	end
    chi2pdftmp = chi2pdf(chi2tmp, mledist.ndata-mledist.npara);
    if(mledist.fiterrors~=1)
        hold on
        plot(chi2tmp, chi2pdftmp, 'r')
        hold off
    end
	xlim([min(chi2tmp) max(chi2tmp)]);
end
title('\chi^2 ')

subplot(3,3,2)
if(mledist.fiterrors==1)
	xbins = linspace(min(mledist.chi2err),max(mledist.chi2err), nbins);
	xbinsdiff = diff(xbins);
	[nc] = histc(mledist.chi2err, xbins);
	nc = nc / sum(nc) / xbinsdiff(1);
	bar(xbins, nc, 'EdgeColor', 'k', 'FaceColor','w');
	arSubplotStyle(gca, labelfontsize, labelfonttype);
	title('\chi_{err}^2')
else
	cla 
	axis off
end

subplot(3,3,3)
if(mledist.fiterrors==1)
	xbins = linspace(min(mledist.chi2fit),max(mledist.chi2fit), nbins);
	xbinsdiff = diff(xbins);
	[nc] = histc(mledist.chi2fit, xbins);
	nc = nc / sum(nc) / xbinsdiff(1);
	bar(xbins, nc, 'EdgeColor', 'k', 'FaceColor','w');
	arSubplotStyle(gca, labelfontsize, labelfonttype);
	title('\chi_{fit}^2 ')
else
	cla 
	axis off
end

g = subplot(3,3,4:9);
try % +ones(100,1)*ar.pTrue(ar.qFit==1)
boxplot(g, mledist.pdiff(:,ar.qFit==1)  , 'orientation', 'horizontal', ...
    'labels', ar.pLabel(ar.qFit==1), 'factordirection', 'list', 'jitter', 1, 'outliersize', 1, ...
    'colors', [0 0 0], 'symbol', 'k.'); 
catch
	boxplot(g, mledist.pdiff(:,ar.qFit==1), 'orientation', 'horizontal', ...
    'labels', ar.pLabel(ar.qFit==1), ...
    'colors', [0 0 0], 'symbol', 'k.');
end
arSubplotStyle(g, labelfontsize, labelfonttype);
title(g, 'difference of truth and estimate')
xlabel(g, 'log_{10}(\theta)-log_{10}(\theta^*)');

if(saveToFile)
    mySaveFigure(h, sprintf('%i_%ix%i', mledist.ntimes, mledist.npoints, mledist.nrep));
    save(sprintf('%s/MLEDist/%i_%ix%i.mat', arSave, ...
        mledist.ntimes, mledist.npoints, mledist.nrep), 'mledist');
end




function h = myRaiseFigure
global mledist
openfigs = get(0,'Children');

figcolor = [1 1 1];

mledist.time = now;

if(isfield(mledist, 'fighandel') && ~isempty(mledist.fighandel) && ...
        mledist.fighandel ~= 0 && ...
        sum(mledist.fighandel==openfigs)>0 && ...
        strcmp(get(mledist.fighandel, 'Name'), 'Distribution of Estimates'))
    
    h = mledist.fighandel;
    figure(h);
else
    h = figure('Name', 'Distribution of Estimates', 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.4 0.4 0.5], 'Color', [1 1 1]);
    set(h,'Color', figcolor);
    mledist.fighandel = h;
end



function mySaveFigure(h, name)
savePath = [arSave '/MLEDist'];

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc2', savePath);
eval(['!ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);



function str = mypath(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');

