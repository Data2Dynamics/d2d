clear all;
close all;
clc;

%% Compile model
arInit

warning off; % switch off warnings -> caution in equilibrium equations
arLoadModel('model_AktPathwayFujita');
warning on;  % switch on warnings

% Load experimental data
for n = 1:6
    arLoadData(['experimentaldata' num2str(n) ],1,'csv');
end

% Compile model
arCompileAll;

% Show condition
arShowDataConditionStructure;

%% Constraint parameters
ar.lb = -8 * ones(size(ar.lb)); % lower parameter bounds
ar.ub =  8 * ones(size(ar.ub)); % upper parameter bounds

%% Numerical settings
ar.config.atol     = 1e-5;
ar.config.rtol     = 1e-5;
ar.config.maxsteps = 1e5;

% Optimizer settings
ar.config.optim.TolFun           = 1e-5;
ar.config.optim.PrecondBandWidth = inf;
ar.config.optim.Display          = 'iter';
ar.config.optim.MaxIter          = 1e4;

%% Optimization
arSimu(true,true,true);
n_fit = 100; % number of starts

if n_fit > 1
    arFitLHS(n_fit); % multistart optimization
else
    arFit;
end

% Print optimization results
arPrint;     % display parameter values
% arPlotChi2s; % waterfall plot

%% Check paper parameters
% arLoadPars('ParamsFujita2010')

%% Visualization of fit
arSimu(false,true,true) % do not compute sensitivities

col_sim = colormap(jet(256));
col_sim = col_sim(1:42:end,:);

figure;
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

for n = 1:6
    plot(ax1,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,3),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax4,ar.model.data(n).tExp,ar.model.data(n).yExp(:,1),ar.model.data(n).yExpStd(:,1),'o--','Color',col_sim(n,:),'LineWidth',1.2);

    plot(ax2,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,4),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax5,ar.model.data(n).tExp,ar.model.data(n).yExp(:,2),ar.model.data(n).yExpStd(:,2),'o--','Color',col_sim(n,:),'LineWidth',1.2);

    plot(ax3,ar.model.condition(n).tFine,ar.model.condition(n).zFineSimu(:,5),'-','Color',col_sim(n,:),'LineWidth',1.5);hold on;
    errorbar(ax6,ar.model.data(n).tExp,ar.model.data(n).yExp(:,3),ar.model.data(n).yExpStd(:,3),'o--','Color',col_sim(n,:),'LineWidth',1.2);
end
time_label = [0:10:60];

xticks(ax1,time_label.*60)
xticklabels(ax1,time_label)
xticks(ax2,time_label.*60)
xticklabels(ax2,time_label)
xticks(ax3,time_label.*60)
xticklabels(ax3,time_label)
xticks(ax4,time_label.*60)
xticklabels(ax4,time_label)
xticks(ax5,time_label.*60)
xticklabels(ax5,time_label)
xticks(ax6,time_label.*60)
xticklabels(ax6,time_label)

ylim(ax1,[0 1.3])
ylim(ax2,[0 1.3])
ylim(ax3,[0 1.3])
ylim(ax4,[0 1.3])
ylim(ax5,[0 1.3])
ylim(ax6,[0 1.3])

xlabel(ax1,'Time [min]')
ylabel(ax1,'Signal [a.u.]')
xlabel(ax2,'Time [min]')
ylabel(ax2,'Signal [a.u.]')
xlabel(ax3,'Time [min]')
ylabel(ax3,'Signal [a.u.]')
xlabel(ax4,'Time [min]')
ylabel(ax4,'Signal [a.u.]')
xlabel(ax5,'Time [min]')
ylabel(ax5,'Signal [a.u.]')
xlabel(ax6,'Time [min]')
ylabel(ax6,'Signal [a.u.]')

title(ax1,'pEGFR - Model')
title(ax2,'pAkt - Model')
title(ax3,'pS6 - Model')
title(ax4,'pEGFR - Data')
title(ax5,'pAkt - Data')
title(ax6,'pS6 - Data')
