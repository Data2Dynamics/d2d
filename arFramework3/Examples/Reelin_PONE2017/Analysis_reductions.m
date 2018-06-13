% Load models & data
%clear all, close all
arInit;
arLoadModel('model_Reelin_beforeReductions');

arLoadData('Reelin_norm', 1);
arLoadData('Reelin_wInh', 1);
arLoadData('Reelin_wInh_Reelin', 1);
arLoadData('Reelin_DR', 1);

arCompileAll;

%Set wide bounds for scaling parameters
arSetPars(ar.pLabel(strncmp('scale',ar.pLabel,5)),ones(1,150)*1,ones(1,150),ones(1,150),ones(1,150)*(-5),ones(1,150)*5);
%Non log bounds for inhibition ratio [0,1]
arSetPars('SFK_inhibition_ratio',-1,1,1,-10,0)

%Structural non-identifiabilities and model reductions
arSetPars({'ApoER2_Dab1_act','init_ApoER2_Dab1','init_pSFK_Int'},[6 0 0],[2 2 2],[1 1 0],[-2 -2 -2],[7 3 3]);

%Use standard deviations from data sheet
ar.config.fiterrors=-1;
ar.config.ploterrors=1;
ar.config.maxsteps = 1e4;

arSetParsPattern('sd_',-1,0,1,-5,3)

arLoadPars('Reelin_nonReduced')

arSave('Reelin_reduction')
ar.config.atol = 1.e-9;
ar.config.rtol = 1.e-9;
arPLEInit
ple([2 6 40])

%PLE of AKT deactivation, dependend on its activation leading to algebraic
%equation relating it to phosphorylated Dab1
plePlot(2)

%PLE of Dab1 trans-phosphorylation open to infinity, linearly dependend on
%SFK inhibition ratio, which leads to introdcution of their product
plePlot(6)

%PLE of potential SFK degradation, open to 0 (-inf on log scale) without 
%dependencies, leading to omission of SFK degradation in model
plePlot(40)


%% Second part of Reelin reductions
arInit
arLoadModel('model_Reelin_withinReductions');

arLoadData('Reelin_norm', 1);
arLoadData('Reelin_wInh', 1);
arLoadData('Reelin_wInh_Reelin', 1);
arLoadData('Reelin_DR', 1);

arCompileAll;

%Set wide bounds for scaling parameters
arSetPars(ar.pLabel(strncmp('scale',ar.pLabel,5)),ones(1,150)*1,ones(1,150),ones(1,150),ones(1,150)*(-5),ones(1,150)*5);
%Non log bounds for inhibition ratio [0,1]
arSetPars('SFK_inhibition_ratio',-1,1,1,-10,0)

%Structural non-identifiabilities and model reductions
arSetPars({'ApoER2_Dab1_act','init_ApoER2_Dab1','init_pSFK_Int'},[6 0 0],[2 2 2],[1 1 0],[-2 -2 -2],[7 3 3]);

%Use standard deviations from data sheet
ar.config.fiterrors=-1;
ar.config.ploterrors=1;
ar.config.maxsteps = 1e4;

arSetParsPattern('sd_',-1,0,1,-5,3)

arLoadPars('Reelin_withinReduced')

arSave('Reelin_within_reduction')
ar.config.atol = 1.e-9;
ar.config.rtol = 1.e-9;
arPLEInit
ple(4)

%PLE of trans-phosph is linearly dependend on
%SFK_release parameter. Summarizing, all three parameters regulating
%binding of Inhibitor to SFKs, its release and its blocking of Dab1
%trans-phosphorylation are correlated. To resolve, SFK_release is expressed
%via trans-phosphorylation and inhibition_ratio in final model
plePlot(4)

%% Third part of Reelin reductions
arInit
arLoadModel('model_Reelin_FINAL');

arLoadData('Reelin_norm', 1);
arLoadData('Reelin_wInh', 1);
arLoadData('Reelin_wInh_Reelin', 1);
arLoadData('Reelin_DR', 1);

arCompileAll;

%Set wide bounds for scaling parameters
arSetPars(ar.pLabel(strncmp('scale',ar.pLabel,5)),ones(1,150)*1,ones(1,150),ones(1,150),ones(1,150)*(-5),ones(1,150)*5);
%Non log bounds for inhibition ratio [0,1]
arSetPars('SFK_inhibition_ratio',-1,1,1,-10,0)

%Structural non-identifiabilities and model reductions
arSetPars({'ApoER2_Dab1_act','init_ApoER2_Dab1','init_pSFK_Int'},[6 0 0],[2 2 2],[1 1 0],[-2 -2 -2],[7 3 3]);

%Use standard deviations from data sheet
ar.config.fiterrors=-1;
ar.config.ploterrors=1;
ar.config.maxsteps = 1e4;

arSetParsPattern('sd_',-1,0,1,-5,3)

arLoadPars('BestFit')
ar.qFit(4) = 1;
arFit
arSave('Reelin_reduction_3')
ar.config.atol = 1.e-9;
ar.config.rtol = 1.e-9;
arPLEInit
ple(4)

%Finally, Dab1 trans-activation parameter is open to infinity, without
%re-optimizations. It is fixed to a high value on the flat end of the
%profile to remain readability of the model
plePlot(4)