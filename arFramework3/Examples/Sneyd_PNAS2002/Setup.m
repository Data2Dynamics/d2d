%% Compile model
arInit
arLoadModel('model_Sneyd_PNAS2002');
arLoadData('dataset_Sneyd_PNAS2002__fig4',1);
arLoadData('dataset_Sneyd_PNAS2002__fig5',1);
arCompileAll;

%% Parameter settings
arSetPars('k1' ,log10(0.64),1,1,-3,5);
arSetPars('k_1',log10(0.04),1,1,-3,5);
arSetPars('k2' ,log10(37.4),1,1,-3,5);
arSetPars('k_2',log10(1.4) ,1,1,-3,5);
arSetPars('k3' ,log10(0.11),1,1,-3,5);
arSetPars('k_3',log10(29.8),1,1,-3,5);
arSetPars('k4' ,log10(4)   ,1,1,-3,5);
arSetPars('k_4',log10(0.54),1,1,-3,5);
arSetPars('l2' ,log10(1.7) ,1,1,-3,5);
arSetPars('l_2',log10(0.8) ,1,1,-3,5);
arSetPars('l4' ,log10(1.7) ,1,1,-3,5);
arSetPars('l_4',log10(2.5) ,1,1,-3,5);
arSetPars('l6' ,log10(4707),1,1,-3,5);
arSetPars('l_6',log10(11.4),1,1,-3,5);

arFitLHS(50)
arPlotChi2s

arPlot
