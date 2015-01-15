% arPlotParameterProfiles(jks)

function arPlotParameterProfiles(jks)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

figure(1); clf;

dchi2 = chi2inv(0.95, 1);

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;

[nrows, ncols] = NtoColsAndRows(length(jks), rowstocols);

count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    mySubplotStyle(g, labelfontsize, labelfonttype)

    plot(ar.ps(:,j), log10(ar.chi2s-min(ar.chi2s)+1), 'xk');
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r--');
    plot(xlim, [0 0], 'k--');
    plot(xlim, log10([dchi2 dchi2]), 'k:');
    hold off
    count = count + 1;
    title(myNameTrafo(ar.pLabel{j}));
    arSpacedAxisLimits;
end



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);


function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);