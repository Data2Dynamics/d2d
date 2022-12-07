%% compile model
arInit

arLoadModel('ABC_model');
arLoadData('ABC_data_BCobs_cov'); % dummy file without data
arCompileAll();

%% set paremeter values
arSetPars('init_A_state',0,1,1,-3,3);
arSetPars('p1',log10(0.05),1,1,-3,3);
arSetPars('p2',log10(0.1),1,1,-3,3);
arSetPars('cov_B',log10(0.3),1,1,-5,log10(0.99));
arSetPars('cov_C',log10(0.3),1,1,-5,log10(0.99));
%arSetPars('cov_B',0,1,0,0,0.99);
%arSetPars('cov_C',0,1,0,0,0.99);
arSetPars('sd_B_au',log10(0.1),1,1,-5,3);
arSetPars('sd_C_au',log10(0.3),1,1,-5,3);

%% simulate unequally spaced data points
nData = 50;
tps1 = linspace(min(ar.model.data.tExp),max(ar.model.data.tExp)/2,ceil(nData*3/4));
tps2 = linspace(max(ar.model.data.tExp)/2,max(ar.model.data.tExp),nData-ceil(nData*3/4));
tps = [tps1,tps2(2:end)];
arSimuData([],[],tps)

ar.p = [0.1067,0.9016,1.1277,-1.4905,0.6118,2.9172,-2.1350];
ar.config.optimizer = 20;
arFit
%ar.config.optimizer = 21;
%arFit
arPLEInit;
ple
