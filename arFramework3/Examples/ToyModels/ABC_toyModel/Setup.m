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
%arSimuData

%Take error 0.1 of data def
ar.config.fiterrors = -1;
arSetPars('sd_B_au',[],2);
arSetPars('sd_C_au',[],2);

