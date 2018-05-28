%% source DDmore (http://repository.ddmore.foundation/model/DDMODEL00000127#Overview)

close all; clc;

Setup_mod127

%%
arFitLHS(50)

arPlotFits
arPlotter
arPrint



%% Profile-Likelihood for population parameters
arPLEInit
ple(1:5)
ple(7:8)
plePlotMulti



%% Profile-Likelihood for the standart deviation of the individual parameters
% local Minimum -> possible to regard this parameter as individual
x = flip((-80:50)/20);
yKgrw1 = pleManual(7,x); % x = flip((-60:60)/40);
yTS0 = pleManual(8,x); % x = flip((-80:80)/40);

plot(x,yKgrw1);
ylabel(string('\chi^2'));
xlabel('log10(pID\_Kgrw1\_sd)');

plot(x,yTS0);
ylabel(string('\chi^2'));
xlabel('log10(pID\_TS0\_sd)');

