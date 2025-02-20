%% Compile model
arInit

warning off; % switch off warnings -> caution in equilibrium equations
arLoadModel('model_ChenMSB2009');
warning on;  % switch on warnings

arLoadData('experimentaldata1',1,'csv');
arLoadData('experimentaldata2',1,'csv');
arLoadData('experimentaldata3',1,'csv');
arLoadData('experimentaldata4',1,'csv');

arCompileAll;

%% Parameter settings
arLoadPars('ParamsChen2009')
ar.config.maxsteps = 5.e4;
ar.qFit(ar.p==0 & ar.qLog10==0) = 2;
arSetPars({'AKT_t','EGFR_t','ERK_t'},[],ones(1,3)*2)
arSetPars('kd115',0,2,0,-5,3)

arFindInputs
arSave('compiled')