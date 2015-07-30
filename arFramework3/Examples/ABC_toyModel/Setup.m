% Load models & data
clear all;
arInit

arLoadModel('ABC_model');
arLoadData('ABC_data_BCobs'); %Data with equidistant observation
%of state B and C for t=0,10,..100
%arLoadData('ABC_data_B_sparseObs'); %Data with sparse observation of
                                 %state B
arCompileAll();

%Simulate new data points
%arSetPars('init_A_star.ate',0,1,1,-3,3);
%arSetPars('p1',log10(0.05),1,1,-3,3);
%arSetPars('p2',log10(0.1),1,1,-3,3);

%Take error 0.1 of data def
ar.config.fiterrors = -1;
arSetPars('sd_B_au',[],2);
arSetPars('sd_C_au',[],2);

arFit();

%Calculate prediction bands for the three states
doPPL(1,1,1:3,linspace(0,100,11),0,1,0.25);

%plot prediction bands
ar.config.ploterrors = -1;
arPlot2