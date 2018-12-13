%% Full model of Boehm et al. Journal of proteome research 2014

close all; clear all; clc;

% initialize full model
Setup

%%
arPlot



% phase plane plot

figure

% experimental data
pADat = ar.model.data.yExp(:,1);
pBDat = ar.model.data.yExp(:,2);
rADat = ar.model.data.yExp(:,3);
ptotDat = (rADat .* pADat + ((100 - rADat) .* pBDat)) / 100;

% full model simulations
pAFull = ar.model.data.yFineSimu(:,1);
pBFull = ar.model.data.yFineSimu(:,2);
rAFull = ar.model.data.yFineSimu(:,3);
ptotFull = (rAFull .* pAFull + ((100 - rAFull) .* pBFull)) / 100;

% plot
plot(ptotDat, rADat, 'ok');
xlabel('STAT5A phosphorylation degree [%]');
ylabel('relative STAT5A amount [%]');
axis([0 100 0 100])
hold on
plot(ptotFull, rAFull, 'k');