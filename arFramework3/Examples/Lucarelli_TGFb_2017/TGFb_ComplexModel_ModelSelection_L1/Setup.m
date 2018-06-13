%%
arInit;

arLoadModel('Complex_ModelSelectionProfile');

% arLoadData('MSPL_mod_complexes_all_hepa16_complex2',1,'xls', 1);
arLoadData('MSPL2_complexes_all_hepa16_no_single_complex',1,'xls', 1);
arLoadData('MSPL2_PL_hepa16_pS2_pS3_time_course',1,'xls', 1);
arLoadData('MSPL2_PL_PN_molecules_per_cell_hepa16',1,'xls', 1);
arLoadData('MSPL2_PL_PN_molecules_per_cell_hepa16_ratios',1,'xls', 1);
arLoadData('PL_PN_TGFbR_hepa16',1,'xls', 1);      
% arLoadData('gene_expression_OE',1,'xls',1);        
arLoadData('MSPL2_complexes_all_hepa16_time_course',1,'xls', 1); 
% % % %arLoadData('gene_expression',[],'xls');       

arCompileAll;


basepath = arSave('compiled');
[pfad,ws_compiled]=fileparts(basepath);
[~,ws_compiled]=fileparts(ar.config.savepath)
% ws_compiled = '20141103T124009_compiled'

ar.lb(strmatch('kdiss_SS',ar.pLabel,'exact')) = -8;

