

function arPlotParameterProfiles

global ar

figure(1)

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;

[nrows, ncols] = NtoColsAndRows(sum(ar.qDynamic==1 & ar.qFit==1), rowstocols);

count = 1;
for j=find(ar.qDynamic==1 & ar.qFit==1)
    g = subplot(nrows,ncols,count);
    mySubplotStyle(g, labelfontsize, labelfonttype)

%     plot(ar.ps(:,j), ar.chi2s, 'xk');
%     ylim([min(ar.chi2s) min(ar.chi2s)+length(ar.p)]);
    semilogy(ar.ps(:,j), ar.chi2s-min(ar.chi2s)+1, 'xk');
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r--');
    hold off
    count = count + 1;
    title(myNameTrafo(ar.pLabel{j}));
end



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);


function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);