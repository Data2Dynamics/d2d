global bdt
load('BDT_allFeatures.mat')

bdt.MODEL_Matrix_bkp2=bdt.MODEL_Matrix;
bdt.RNA_Matrix_bkp2 = bdt.RNA_Matrix;

% w/o mutational info
bdt.RNA_Matrix(:,6:7) = NaN;
bdt.MODEL_Matrix(:,11:30) = NaN;

doBDT_calculations(0.3,1:2)

bdt.RNA_Matrix = bdt.RNA_Matrix_bkp2;
bdt.MODEL_Matrix = bdt.MODEL_Matrix_bkp2;

RNA_noMut = plot_Conf('RNA_mut_wHGF_0_3',0,'ALL')
MODEL_noMut = plot_Conf('MODEL_mut_wHGF_0_3',0,'ALL')

save('BDT_no_mutation.mat','bdt','bdt_figures')

% MODEL: only SS
bdt.MODEL_Matrix(:,2:10) = NaN;
bdt.MODEL_Matrix(:,1) = 1;
bdt.MODEL_Matrix(:,[15:18 24:28]) = NaN;

doBDT_calculations(0.3,1)
bdt.MODEL_Matrix = bdt.MODEL_Matrix_bkp2;

MODEL_SS = plot_Conf('MODEL_mut_wHGF_0_3',0,'ALL')
save('BDT_quasiSS.mat','bdt','bdt_figures')

% MODEL: only FoldC
bdt.MODEL_Matrix(:,2:10) = NaN;
bdt.MODEL_Matrix(:,1) = 1;
bdt.MODEL_Matrix(:,[11:14 19:23]) = NaN;

doBDT_calculations(0.3,1)
bdt.MODEL_Matrix = bdt.MODEL_Matrix_bkp2;

MODEL_FoldC = plot_Conf('MODEL_mut_wHGF_0_3',0,'ALL')
save('BDT_FoldC.mat','bdt','bdt_figures')

% MODEL: All
doBDT_calculations(0.3,1)
bdt.MODEL_Matrix = bdt.MODEL_Matrix_bkp2;

MODEL_ALL = plot_Conf('MODEL_mut_wHGF_0_3',0,'ALL')
save('BDT_allFeatures.mat','bdt','bdt_figures')
