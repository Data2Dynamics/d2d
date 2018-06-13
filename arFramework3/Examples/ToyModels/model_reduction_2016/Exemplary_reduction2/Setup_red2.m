% Load models & data
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
 
