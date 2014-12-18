for i=1:50:length(ar.p)
    figure;
    if i+50 > length(ar.p)
        xdata = linspace(i,length(ar.p),length(ar.p)-i+1);
        errorbar(xdata,ar.mean(i:length(ar.p)),ar.std(i:length(ar.p)),'b.');
        xlim([i length(ar.p)]);
        ylim([-5.5 3.5]);
        set(gca,'tickdir','out');
        hold on
        plot(xdata,ar.p(i:length(ar.p)),'r*');
        set(gca,'XTick',xdata);
        set(gca,'XTickLabel',ar.pLabel(i:length(ar.p)));
        view(gca,[90 90]);
        hold off
    else
        xdata = linspace(i,i+49,50);
        errorbar(xdata,ar.mean(i:i+49),ar.std(i:i+49),'b.');
        xlim([i i+49]);
        ylim([-5.5 3.5]);
        set(gca,'tickdir','out');
        hold on
        plot(xdata,ar.p(i:i+49),'r*');
        set(gca,'XTick',xdata);
        set(gca,'XTickLabel',ar.pLabel(i:i+49));
        view(gca,[90 90]);
        hold off
    end
end


% figure;
% hold on
% errorbar(linspace(1,length(ar.p),length(ar.p)),ar.mean,ar.std,'b.');
% xlim([0 length(ar.p)]);
% plot(linspace(1,length(ar.p),length(ar.p)),ar.p,'r*');
% hold off