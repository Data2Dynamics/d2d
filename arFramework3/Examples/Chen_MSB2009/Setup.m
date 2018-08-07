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
arFindInputs

%% Parameter settings
arLoadPars('ParamsChen2009')
ar.config.maxsteps = 5.e4;
