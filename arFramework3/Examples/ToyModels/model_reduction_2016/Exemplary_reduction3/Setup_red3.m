% Load models & data
clear all;
arInit

%Load full or reduced model
arLoadModel('model_expl3_full');
% arLoadModel('model_expl3_red');

%Load data set
arLoadData('data_expl3');
arCompileAll();
%Calibrate model, find global optimum
arFitLHS(10)
%Print parameters
arPrint

%Save model
arSave

%for the full model, calculate PLE of k_d,Z
arPLEInit
ple(8,200,1.e-2)

%for reduced model, calculate PLE of scaling factor between Y and Z
% ple(7)