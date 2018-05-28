%% source: DDmore http://repository.ddmore.foundation/model/DDMODEL00000103#Overview
% data103_IDall: data file from DDmore, 569 individuals, a lot with <5 mesurements,
%   takes very long to compute, only one individual parameter possible
% data103_IDlong: data of all individuals with at least 5 mesurements
% data103_IDshort: data of idividual 508-523, similar results to the long file,
%   significant shorter computation time

close all; clc;

Setup_mod103

%%

arFitLHS(100)

arPlotFits
arPlotter
arPrint(1:8) % only population parameters


%% Profile-Likelihood of the population parameters
arPLEInit
ple(1:8)
plePlotMulti


%% Profile-Likelihood for the standard deviation of the individual parameters
% KDE, Kout, Slop are declared as individual on DDmore,
% if there is a local minimum in their Profile-Likelihood, set the prior
% that the fit converges into that local minimum
% set the prior for one parameter then do another LHS fit, set the next
% prior
x = flip((-60:60)/40);
yKDEr = pleManual(6,x);
yKDEl = pleManual(6,flip(x));
yKDE = min(yKDEr,flip(yKDEl));
plot(x,yKDE);
ylabel(string('\chi^2_{PL}'));
xlabel('log_{10}(KDE_{sd})');

x = flip((-80:80)/40);
yKoutr = pleManual(7,x);
yKoutl = pleManual(7,flip(x));
yKout = min(yKoutr,flip(yKoutl));
plot(x,yKout);
ylabel(string('\chi^2_{PL}'));
xlabel('log_{10}(Kout_{sd})');

x = flip((-40:40)/40);
ySlopr = pleManual(8,x);
ySlopl = pleManual(8,flip(x));
ySlop = min(ySlopr,flip(ySlopl));
plot(x,ySlop);
ylabel(string('\chi^2_{PL}'));
xlabel('log_{10}(Slop_{sd})');



