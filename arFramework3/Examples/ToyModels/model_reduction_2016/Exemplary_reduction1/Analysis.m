
Setup_red1

%%Calibrate the model
arFit;
arPrint
%Save the model
arSave

%For the full model, calculate the profile of k1
arPLEInit
ple(1,150)