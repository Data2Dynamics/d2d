% plot mcmc chains

function arPlotMCMCParVsChi(jks, Nthinning, popt, chi2opt)

global ar

if(~exist('Nthinning','var'))
    Nthinning = 1;
end
ps = ar.ps;
chi2s = ar.chi2s;
if(Nthinning>1)
    ps = ps(mod(1:size(ps,1),Nthinning)==1,:);
    chi2s = chi2s(mod(1:length(chi2s),Nthinning)==1,:);
end

labelfontsize = 12;

h = myRaiseFigure('mcmc chains');
set(h, 'Color', [1 1 1]);

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end

rowstocols = 0.5; %0.7; 0.45;

jks = jks(ar.qFit(jks)==1);

[nrows, ncols] = NtoColsAndRows(length(jks), rowstocols);

count = 1;
for jk=jks
    
    xlimtmp2 = (max(ps(:,jk))-min(ps(:,jk)))*0.05;
    if(xlimtmp2>0)
        xlimtmp = [min(ps(:,jk))-xlimtmp2 max(ps(:,jk))+xlimtmp2];
    end
    
    g = subplot(nrows, ncols, count);
    set(g, 'FontSize', labelfontsize);
    set(g, 'FontName', 'TimesNewRoman');

    plot(ps(:,jk), chi2s, 'kx', 'MarkerSize', 1)
    hold on
    colors = jet(length(popt));
    for j=1:length(popt)
        plot(popt{j}(jk), chi2opt(j), '*', 'Color', colors(j,:));
    end
    hold off

    
    xlim([xlimtmp(1)-xlimtmp2*0.05 xlimtmp(2)+xlimtmp2*0.05]);
    %ylim([1 size(ps,1)]);
    title(myNameTrafo(ar.pLabel{jk}))
    
    count = count + 1;
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

function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);

