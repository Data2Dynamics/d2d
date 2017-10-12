%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   Setup_BDT_calc                          %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Perform calculation of BDTs and efficiency plots
%  Also, TCGA predictions and correlations to ligands can be calculated

%Load BDT matrices directly or from model folder?
load('BDT_plus_TCGA_data.mat')
%Load from model folder (calculate in PanRTK_final_forBDT folder)
%load('../PanRTK_final_forBDT/BDT_matrices_new.mat')
global bdt
global bdt_figures

run_regression = 0;
%Run calculations with n% of cell lines used for testing on both model and
%RTK expression data sets
doBDT_calculations(0.3, 1:2)

%get efficiency for all ligands and their combination for model and RNA,
%produce plots
plot_BDT_einzel_RNAvsModel(1)


% %Comment out to train BDT on all cell lines for TCGA output
% %Load own data set (calculate in PanRTK_final_forBDT folder)
%load('../PanRTK_final_forBDT/Own_BDT_plus_TCGA_matrices.mat')

addpath('distributionPlot')
Mat_MODEL_mut = bdt.MODEL_Matrix;
Mat_MODEL_mut(isnan(Mat_MODEL_mut(:,end)),:)=[];
Mat_MODEL_mut((Mat_MODEL_mut(:,bdt.nr_model+3)==1),:) = NaN;

a = TreeBagger(500,Mat_MODEL_mut(:,1:bdt.nr_model),Mat_MODEL_mut(:,end),'SampleWithReplacement','on','Method','classification','oobvarimp','on','MinLeaf',10,'FBoot',0.9);

%Get model predictions on TCGA patient data
plot_TCGA_predictions

%Get correlation to ligand expression
plot_TCGA_corr_Ligand(run_regression)