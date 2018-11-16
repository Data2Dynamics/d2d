% Load models & data
%clear all, close all
arInit;
%State if the initial or final model is to be loaded
load_final = 1;

if(load_final)
    arLoadModel('model_Reelin_FINAL');
else
    arLoadModel('model_initial_Reelin');
end

arLoadData('Reelin_norm', 1);
arLoadData('Reelin_wInh', 1);
arLoadData('Reelin_wInh_Reelin', 1);
arLoadData('Reelin_DR', 1);
arLoadData('Reelin_EC50', 1);

arCompileAll;
arFindInputs;
%Set wide bounds for scaling parameters
arSetPars(ar.pLabel(strncmp('scale',ar.pLabel,5)),ones(1,150)*1,ones(1,150),ones(1,150),ones(1,150)*(-5),ones(1,150)*5);
%Non log bounds for inhibition ratio [0,1]
arSetPars('SFK_inhibition_ratio',-1,1,1,-10,0)

%Structural non-identifiabilities and model reductions
if(load_final)
    arSetPars({'ApoER2_Dab1_act','init_ApoER2_Dab1'},[6 0],[2 2],ones(1,2),[-2 -2],[7 3]);
else
    arSetPars({'init_pbApoER2','init_pSFK_Int','offset_Dab1_ExpReelinStim','offset_pDab1_ExpReelin_DR','SFK_activation'},[0  0 0 0 3],[2 2 2 2 2],[1 1 0 0 1],[-3 -3 -3 -3 -3],[3 3 3 3 5])
end

%Use standard deviations from data sheet
ar.config.fiterrors=-1;
ar.config.ploterrors=1;
ar.config.maxsteps = 1e4;

arSetParsPattern('sd_',-1,0,1,-5,3)

if(load_final)
    arLoadPars('BestFit')
else
    arLoadPars('Initial_BestFit')
end

%Fix observational parameters without corresponding data
arSetPars({'offset_tApo_ExpEC50R','offset_tApo_ExpReelin_DR','offset_tApo_ExpReelin_Inh','offset_tApo_ExponlyInh','offset_tSFK_ExpEC50R','offset_tSFK_ExpReelin_DR','offset_tSFK_ExpReelin_Inh','offset_tSFK_ExponlyInh','scale_tApo__ExpEC50R','scale_tApo__ExpReelin_DR','scale_tApo__ExpReelin_Inh','scale_tApo__ExponlyInh','scale_tSFK_ExpEC50R','scale_tSFK_ExpReelin_DR','scale_tSFK_ExpReelin_Inh','scale_tSFK_ExponlyInh',},zeros(1,16),ones(1,16)*2,zeros(1,16),ones(1,16)*-5,ones(1,16)*3)
arSetPars('offset_tSFK_ExpReelinStim',[],2,1,-3,3)

%No log plotting
ar.model(1).qPlotYs(:) = 1;
for i = 1:length(ar.model(1).data)
    ar.model(1).data(i).logplotting(:) = 0;
end

% %Part for prediction profile calculations
% if(load_final)
%     arLoadPars('BestFit_VPLPars')
% else
%     arLoadPars('Initial_BestFit_VPLPars')
% end
% %store and del EC50 data
% yStd_bkp = ar.model(1).data(18).yExpStd;
% yExp_bkp = ar.model(1).data(18).yExp;
% ar.model(1).data(18).yExpStd(:) = NaN;
% ar.model(1).data(18).yExp(:) = NaN;
% 
% %do predictions
% PPL_options('alpha_level',0.34,'onlyProfile',false);
% doPPL(1,18,[1:4],[1 5 10 15 25 45 60 100 150 240],1)
% 
% ar.model(1).data(18).yExpStd(:) = yStd_bkp;
% ar.model(1).data(18).yExp(:) = yExp_bkp;
% ar.model(1).qPlotYs(1:end-1) = 0;
% arPlot2
