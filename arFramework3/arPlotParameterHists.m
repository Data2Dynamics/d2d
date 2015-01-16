

function arPlotParameterHists(ps, jks, nbins)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

if(~exist('nbins','var'))
    nbins = 50;
end

figure(1); clf;

[nrows, ncols] = arNtoColsAndRows(length(jks));

count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    arSubplotStyle(g)
    hist(ps(:,j), nbins);
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r');
    hold off
    count = count + 1;
    title(arNameTrafo(ar.pLabel{j}));
end







