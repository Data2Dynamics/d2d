% Plot stochastic noise

function arPlotSSAnoise(m, c, scaling)

global ar;

figure(1234);
colors = lines(size(ar.model(m).condition(c).xSSA,2));
for jx = 1:size(ar.model(m).condition(c).xSSA,2)
    qt = ar.model(m).condition(c).tFine > ar.model(m).condition(c).tFine(end)*0.05;
    x = ar.model(m).condition(c).xFineSimu(qt,jx)*scaling(jx);
    
    xmean = mean(squeeze(ar.model(m).condition(c).xSSA(qt,jx,:)*scaling(jx)), 2);
    xstd = std(squeeze(ar.model(m).condition(c).xSSA(qt,jx,:)*scaling(jx)), 1, 2);
    
    [xmean, xi] = sort(xmean);
    xstd = xstd(xi,:);
    
    qnonzero = xmean~=0 & xstd>1e-6;
    xmean = xmean(qnonzero);
    xstd = xstd(qnonzero,:);
    
    loglog(xmean, xstd./xmean, 'Color', colors(jx,:));
    hold on
end

hold off

legend(ar.model(m).x)
title(sprintf('amount of stochastic noise (%i SSA runs)', ar.model(m).condition(c).nSSA))
ylabel('standard deviation [%]')
xlabel('mean signal [number of molecules]')
