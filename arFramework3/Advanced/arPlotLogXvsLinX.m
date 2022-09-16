function arPlotLogXvsLinX

global ar

% arLoadLatest('BothScales')
% if length(ar.model)<2
%     return 
% end
% arPlot
% close all

figure % new figure to not plot into d2d Subplots

[~,ia,ib] = intersect(ar.model(1).x,ar.model(2).z,'stable');
maxdiff = [];
maxreldiff = [];
for i=1:length(ia)
    for j = 1:length(ar.model(1).condition)
        tFine = ar.model(1).condition(j).tFine;
        xlog = ar.model(1).condition(j).xFineSimu(:,ia(i));
        zlog = interp1(ar.model(2).condition(j).tFine,ar.model(2).condition(j).zFineSimu(:,ib(i)),tFine);
        loglog(xlog,...
            zlog,...
            '.','MarkerSize',12)
        hold on
        maxdiff(i) = max(abs(xlog-zlog));
        maxreldiff(i) = max(abs(xlog-zlog)./((xlog+zlog)/2));
        set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')+1);
    end
end
maxdiff
maxreldiff
myxlim = xlim;
plot(myxlim,myxlim,'k')
axis tight
xlabel('standard x')
ylabel('log(x)');
legend(ar.model(1).x{:},'Location','SouthEast','Interpreter','none');
C = strsplit(ar.info.path,'/');
title(strrep(C{end},'_','\_'));
saveas(gcf,'Comparison_LogandLinStates.png')