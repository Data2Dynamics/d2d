%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   comp_boxplot_validation                 %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Create boxplot comparing random receptor/surcace levels with intial fit
%  of signaling model to calibration/validation cell lines
%
%

nr_celllines = ar.val.plots+ar.val.nr_val;

chi2_cal = nansum(ar.val.chi2(1:ar.val.plots))/nansum(ar.val.ndata(1:ar.val.plots));
chi2_val = nansum(ar.val.chi2(ar.val.plots+1:nr_celllines))/nansum(ar.val.ndata(ar.val.plots+1:nr_celllines));
nr_rnd = floor(length(ar.val.chi2(ar.val.plots:end))/nr_celllines);

for i=2:nr_rnd
    chi2_valrnd(i-1) = nansum(ar.val.chi2((i-1)*nr_celllines+ar.val.plots+1:(i-1)*nr_celllines+nr_celllines))/nansum(ar.val.ndata((i-1)*nr_celllines+ar.val.plots+1:(i-1)*nr_celllines+nr_celllines));
    chi2_rnd(i-1) = nansum(ar.val.chi2((i-1)*nr_celllines+1:(i-1)*nr_celllines+ar.val.plots))/nansum(ar.val.ndata((i-1)*nr_celllines+1:(i-1)*nr_celllines+ar.val.plots));
end
while(chi2_rnd(1)==0)
    chi2_rnd(1) = [];
end
while(chi2_rnd(1)==0)
    chi2_valrnd(1) = [];
end
chi2_valrnd([1:3 12])=NaN;
% boxplot([chi2_rnd' chi2_valrnd'],[repmat('Calibration',length(chi2_rnd),1); repmat('Validation ',length(chi2_valrnd),1)])%,'boxstyle','filled','plotstyle','compact'
boxplot(chi2_rnd',repmat('Calibration',length(chi2_rnd),1));%,'boxstyle','filled','plotstyle','compact'

hold on
set(gca,'YLim',[2.8 5.5])
plot(1,[chi2_cal],'.','Color',[0.55 0 0],'MarkerSize',20,'MarkerFaceColor',[0.55 0 0])
plot(1,[chi2_val],'.','Color',[0 0 0.55],'MarkerSize',20,'MarkerFaceColor',[0.55 0 0])

[~,p_val] = ttest(chi2_rnd,chi2_cal);
[~,p_val2] = ttest(chi2_valrnd(~isnan(chi2_valrnd)),chi2_val);

ylabel('\chi^2 /N','Interpreter','Tex');
title('Goodness of fit for calibration and validation cell lines');
set(gcf,'Color','w')
text(0.7,5,['p-value = ' num2str(p_val)])
text(1.7,5,['p-value = ' num2str(p_val2)])

mean(chi2_rnd)
std(chi2_rnd)

mean(chi2_valrnd)
std(chi2_valrnd)