arInit

arLoadModel('ABC_model');
arLoadData('ABC_data_2'); 
arCompileAll();

% Fix Errors from data sheet
ar.config.fiterrors = -1;
arSetPars('sd_A_au',[],2);
arSetPars('sd_B_au',[],2);
arSetPars('sd_C_au',[],2);

% % True parameters: 
% arSetPars('init_A_state',0,1,1,-3,3);
% arSetPars('p1',log10(0.05),1,1,-3,3);
% arSetPars('p2',log10(0.1),1,1,-3,3);

arLoadPars('BestFit');