%% source: DDmore (http://repository.ddmore.foundation/model/DDMODEL00000118#Overview)
close all; clc;

Setup_mod118

%%
arFitLHS(500)

arPlotFits
arPlotter
arPrint(1:12) % population parameters only


%% Profile-Likelihood for population parameters
% global minimum (d.h. V2 = 0) in conflict with three compartment assumption
arPLEInit
ple(1:7)
ar.lb(9) = -0.9;
ple(9)
plePlotMulti


%% Profile-Likelihood for the standart deviation of the individual parameters
% local Minimum -> possible to regard this parameter as individual
x = flip((-60:40)/40);
yCL = pleManual(9,x);
yQ = pleManual(11,x);
yV2 = pleManual(12,x);

plot(x,yCL);
ylabel(string('\chi^2'));
xlabel('log10(pID\_CL\_sd)');

plot(x,yQ);
ylabel(string('\chi^2'));
xlabel('log10(pID\_Q\_sd)');

plot(x,yV2);
ylabel(string('\chi^2'));
xlabel('log10(pID\_V2\_sd)');



