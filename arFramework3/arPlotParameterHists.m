

function arPlotParameterHists(ps, jks, nbins)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

if(~exist('nbins','var'))
    nbins = 50;
end

figure(1); clf;

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;

[nrows, ncols] = NtoColsAndRows(length(jks), rowstocols);

count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    arSubplotStyle(g, labelfontsize, labelfonttype)
    hist(ps(:,j), nbins);
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r');
    hold off
    count = count + 1;
    title(myNameTrafo(ar.pLabel{j}));
end



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);


function str = myNameTrafo(str)
str = strrep(str, '_', '\_');


