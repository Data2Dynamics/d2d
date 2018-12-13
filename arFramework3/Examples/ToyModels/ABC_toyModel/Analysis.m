% Load models & data
clear all;

Setup

arFit();

%Calculate prediction bands for a state
PPL_options('Integrate',true,'Stepsize',0.5)
arPPL(1,1,1,linspace(0,100,11),1);

%plot prediction bands
arPlot2
