% Load models & data
Setup_red2

%%
%Calibrate the model
arFit
%Print the parameters
arPrint
%Save the model
arSave

%calculate profile of k_4 for full model or k_5 for reduced model
arPLEInit
ple(6)