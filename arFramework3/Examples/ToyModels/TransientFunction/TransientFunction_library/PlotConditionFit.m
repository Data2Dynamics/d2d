% dat

function PlotConditionFit(dat)
global ar

if length(ar.model)>1
    error('Use a model where only the transient fucntion is available.')
end
if length(ar.model.data)>1
    error('Use a model where only a single data set is available.')
end

% close all
clf('reset')

tmin = min(ar.model.data.tExp(~isnan(ar.model.data.yExp)));
tmax = max(ar.model.data.tExp(~isnan(ar.model.data.yExp)));

ind = find(ar.model.data.tFine<=tmax & ar.model.data.tFine>=tmin);
patch([ar.model.data.tFine(ind);ar.model.data.tFine(ind(end:-1:1))],...
    [ar.model.data.yFineSimu(ind)+ar.model.data.ystdFineSimu(ind);ar.model.data.yFineSimu(ind(end:-1:1))-ar.model.data.ystdFineSimu(ind(end:-1:1))],...
    zeros(length(ind)*2,1)-2*eps,'FaceColor',0.7*ones(1,3),'EdgeColor',0.7*ones(1,3));
hold on
h1=plot(ar.model.data.tFine(ind),ar.model.data.yFineSimu(ind),'k','LineWidth',2);
hold on


ind = dat.tExp<=tmax & dat.tExp>=tmin;
h2=plot(dat.tExp(ind),dat.yExp,'b','LineWidth',2);

set(gca,'XLim',[tmin,tmax]+[-1,1]*0.02*(tmax-tmin));
legend([h1,h2],'TF','ODE')

