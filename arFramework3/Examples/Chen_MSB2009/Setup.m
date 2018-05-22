clear all;
close all;
clc;

%% Compile model
arInit

warning off; % switch off warnings -> caution in equilibrium equations
arLoadModel('model_ChenMSB2009');
warning on;  % switch on warnings

arLoadData('experimentaldata1',1,'csv');
arLoadData('experimentaldata2',1,'csv');
arLoadData('experimentaldata3',1,'csv');
arLoadData('experimentaldata4',1,'csv');

arCompileAll;

arShowDataConditionStructure;

%% Parameter settings
load all.mat
arLoadPars('ParamsChen2009')

%% Visualization of fit
arSimu(false,true,true) % do not compute sensitivities

col_sim = colormap(hsv(4));

ax1 = subplot(2,3,1); hold all;
axis square
ax2 = subplot(2,3,2); hold all;
axis square
ax3 = subplot(2,3,3); hold all;
axis square
ax4 = subplot(2,3,4); hold all;
axis square
ax5 = subplot(2,3,5); hold all;
axis square
ax6 = subplot(2,3,6); hold all;
axis square

for n = 1:4
    plot(ax1,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,1),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax4,ar.model.data(n).tExp,ar.model.data(n).yExp(:,1),ar.model.data(n).yExpStd(:,1),'o--','Color',col_sim(n,:),'LineWidth',1.2);

    plot(ax2,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,2),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax5,ar.model.data(n).tExp,ar.model.data(n).yExp(:,2),ar.model.data(n).yExpStd(:,2),'o--','Color',col_sim(n,:),'LineWidth',1.2);

    plot(ax3,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,3),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax6,ar.model.data(n).tExp,ar.model.data(n).yExp(:,3),ar.model.data(n).yExpStd(:,3),'o--','Color',col_sim(n,:),'LineWidth',1.2);
end
time_label = [0:10:60];

ylim(ax1,[0 1.2])
ylim(ax2,[0 1.2])
ylim(ax3,[0 1.2])

ylim(ax4,[0 120])
ylim(ax5,[0 120])
ylim(ax6,[0 120])

xlabel(ax1,'Time [a.u.]')
ylabel(ax1,'Signal [a.u.]')
xlabel(ax2,'Time [a.u.]')
ylabel(ax2,'Signal [a.u.]')
xlabel(ax3,'Time [a.u.]')
ylabel(ax3,'Signal [a.u.]')
xlabel(ax4,'Time [a.u.]')
ylabel(ax4,'Signal [a.u.]')
xlabel(ax5,'Time [a.u.]')
ylabel(ax5,'Signal [a.u.]')
xlabel(ax6,'Time [a.u.]')
ylabel(ax6,'Signal [a.u.]')

title(ax1,'pErbB1 - Model')
title(ax2,'pErk - Model')
title(ax3,'pAkt - Model')
title(ax4,'pErbB1 - Data')
title(ax5,'pErk - Data')
title(ax6,'pAkt - Data')
