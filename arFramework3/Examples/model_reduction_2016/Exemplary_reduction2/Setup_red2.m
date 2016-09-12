% Load models & data
clear all;
arInit
%Load either the full or the reduced model
arLoadModel('model_expl2_full');
% arLoadModel('model_expl2_red');

%Load the data set
arLoadData('data_expl2');
arCompileAll();

arSetPars('init_pX_state',0,2,0,-5,3);
arSetPars('init_ppX_state',0,2,0,-5,3);
 ar.p(5) = -3;
 ar.qFit(5) = 2;
 
%Calibrate the model
arFit
%Print the parameters
arPrint
%Save the model
arSave

%calculate profile of k_4 for full model or k_5 for reduced model
arPLEInit
ple(6)