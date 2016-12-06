%%
arInit;
arLoadModel('model_template');
arLoadData('data_b10');
arLoadData('data_betacar');
arLoadData('data_betacar_zea');
arLoadData('data_betacry');
arLoadData('data_ohb10');
arLoadData('data_zea');
arCompileAll;

%%
ar.config.ploterrors=1;
ar.config.fiterrors=-1;

