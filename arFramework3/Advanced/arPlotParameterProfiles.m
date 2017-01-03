% arPlotParameterProfiles(jks)

function arPlotParameterProfiles(jks)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

figure(1); clf;

dchi2 = chi2inv(0.95, 1);

[nrows, ncols] = arNtoColsAndRows(length(jks));

count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    arSubplotStyle(g)

    plot(ar.ps(:,j), log10(ar.chi2s-min(ar.chi2s)+1), 'xk');
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r--');
    plot(xlim, [0 0], 'k--');
    plot(xlim, log10([dchi2 dchi2]), 'k:');
    hold off
    count = count + 1;
    title(arNameTrafo(ar.pLabel{j}));
    arSpacedAxisLimits;
end


