arInit;

% Config settings:
% ar.config.optim.TolFun = 1e-9;
% ar.config.optim.TolX = 0; 

arLoadModel('Complex_ModelSelectionProfile'); 

arLoadData('MSPL2_complexes_all_hepa16_no_single_complex',1,'xls', 1);
arLoadData('MSPL2_PL_hepa16_pS2_pS3_time_course',1,'xls', 1);
arLoadData('MSPL2_PL_PN_molecules_per_cell_hepa16',1,'xls', 1);
arLoadData('MSPL2_PL_PN_molecules_per_cell_hepa16_ratios',1,'xls', 1);
arLoadData('PL_PN_TGFbR_hepa16',1,'xls', 1);      
arLoadData('gene_expression_OE',1,'xls',1);
arLoadData('MSPL2_complexes_all_hepa16_time_course',1,'xls', 1); 

global ar
ar.config.no_optimization = 1;
arCompileAll;

arLoadPars('BestFit');

arPrint
