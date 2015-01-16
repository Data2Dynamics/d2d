% check distribution of MLE
%
% arMLEDistCompare(mledist, mlelabels, saveToFile)
%
% plotChi2dist      [true]
% saveToFile        [false]

function arMLEDistCompare(mledist, mlelabels, saveToFile)

global ar

if(~exist('saveToFile','var'))
    saveToFile = false;
end

% constants
labelfontsize = 8;
labelfonttype = 'TimesNewRoman';

hfig = myRaiseFigure;

for j=1:length(mledist)
    datachi2(:,j) = mledist{j}.chi2; %#ok<*AGROW>
    datachi2err(:,j) = mledist{j}.chi2err;
    datachi2fit(:,j) = mledist{j}.chi2fit;
end

colors = lines(length(mledist));

datapdiff = [];
datalabels = {};
datapos = [];
pcount = 1;
morecolors = [];
for j=find(ar.qFit==1)
    for jj=1:length(mledist)
        datapos(end+1) = jj + (pcount-1)*(length(mledist)+2);
        morecolors(end+1,:) = colors(jj,:);
        
        datapdiff(:,end+1) = mledist{jj}.pdiff(:,j);
        if(jj == ceil(length(mledist)/2))
            datalabels(end+1) = ar.pLabel(j);
        else
            datalabels{end+1} = '';
        end
    end
    pcount = pcount + 1;
end

nbins = 50;

for j=1:length(mledist)
    subplot(2+length(mledist),3,(j-1)*3+1)
    
    xbins = linspace(min(datachi2(:,j)),max(datachi2(:,j)), nbins);
    xbinsdiff = diff(xbins);
    
    [nc] = histc(datachi2(:,j), xbins);
    nc = nc / sum(nc) / xbinsdiff(1);
    bar(xbins, nc, 'EdgeColor', colors(j,:), 'FaceColor','w');
    hold on
    
    chi2tmp = linspace(min(datachi2(:,j)),max(datachi2(:,j)), 100);
    chi2pdftmp(:,j) = chi2pdf(chi2tmp, mledist{j}.ndata-mledist{j}.npara);
    plot(chi2tmp, chi2pdftmp(:,j), 'Color', colors(j,:))
    hold off
    
    arSubplotStyle(gca, labelfontsize, labelfonttype);
    
    if(j==1)
        title('\chi^2 ')
    end
    ylabel(mlelabels{j});
end

for j=1:length(mledist)
    subplot(2+length(mledist),3,(j-1)*3+2)
    
    xbins = linspace(min(datachi2err(:,j)),max(datachi2err(:,j)), nbins);
    xbinsdiff = diff(xbins);
    
    nc = histc(datachi2err(:,j), xbins);
    nc = nc / sum(nc) / xbinsdiff(1);
    bar(xbins, nc, 'EdgeColor', colors(j,:), 'FaceColor','w');
    
    arSubplotStyle(gca, labelfontsize, labelfonttype);
    
    if(j==1)
        title('\chi_{err}^2 ')
    end
end

for j=1:length(mledist)
    subplot(2+length(mledist),3,(j-1)*3+3)
    
    xbins = linspace(min(datachi2err(:,j)),max(datachi2fit(:,j)), nbins);
    xbinsdiff = diff(xbins);
    
    nc = histc(datachi2err(:,j), xbins);
    nc = nc / sum(nc) / xbinsdiff(1);
    bar(xbins, nc, 'EdgeColor', colors(j,:), 'FaceColor','w');
    
    arSubplotStyle(gca, labelfontsize, labelfonttype);
    
    if(j==1)
        title('\chi_{fit}^2 ')
    end
end

g = subplot(2+length(mledist),3,(length(mledist)*3+1):((2+length(mledist))*3));
boxplot(g, datapdiff, 'orientation', 'horizontal', ...
    'labels', datalabels, 'outliersize', 1, 'jitter', 1, 'factordirection', 'list', ...
    'colors', morecolors, 'symbol', 'k.', ...
    'positions', datapos);
arSubplotStyle(g, labelfontsize, labelfonttype);
title(g, 'difference of truth and estimate')
xlabel(g, 'log_{10}(\theta-\theta^*)');

if(saveToFile)
    mySaveFigure(hfig, 'mledist_compare');
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
savePath = [arSave '/MLECompare'];

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

