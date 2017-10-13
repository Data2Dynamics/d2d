% Load models & data
%clear all, close all
arInit;

arLoadModel('model_Reelin_FINAL');

arLoadData('Reelin_Corrected_1210_Dab1Src_Log', 1);
arLoadData('Reelin_Corrected_2303_wInh_Log', 1);
arLoadData('Reelin_Corrected_1210_wInhR_Log', 1);
arLoadData('Reelin_Corrected_1410_DR_Log', 1);
arLoadData('Reelin_Corrected_1412_EC50_log', 1);

arCompileAll;

%Set wide bounds for scaling parameters
arSetPars(ar.pLabel(strncmp('scale',ar.pLabel,5)),ones(1,150)*1,ones(1,150),ones(1,150),ones(1,150)*(-5),ones(1,150)*5);
%Non log bounds for inhibition ratio [0,1]
arSetPars('SFK_inhibition_ratio',-1,1,1,-10,0)

%Structural non-identifiabilities and model reductions
arSetPars({'ApoER2_Dab1_act','init_ApoER2_Dab1'},[6 0],[2 2],ones(1,2),[-2 -2],[7 3]);

%Use standard deviations from data sheet
ar.config.fiterrors=-1;
ar.config.ploterrors=1;
ar.config.maxsteps = 1e4;

arSetParsPattern('sd_',-1,0,1,-5,3)

arLoadPars('BestFit')

%No log plotting
ar.model(1).qPlotYs(:) = 1;
for i = 1:length(ar.model(1).data)
    ar.model(1).data(i).logplotting(:) = 0;
end