
arInit;

arLoadModel('combined_model');

%% Incucyte Data

arLoadData('Incucyte/U2OS_Gem', 'combined_model', [], true);
arLoadData('Incucyte/U2OS_Dox', 'combined_model', [], true);
arLoadData('Incucyte/U2OS_SN38', 'combined_model', [], true);

%% WB data (time course)
arLoadData('WB/TC_Dox', 'combined_model', [], true);
arLoadData('WB/TC_SN38', 'combined_model', [], true);
arLoadData('WB/TC_Gem', 'combined_model', [], true);

%% WB data (6h dose response)
arLoadData('WB/DR_Dox', 'combined_model', [], true);
arLoadData('WB/DR_SN38', 'combined_model', [], true);
arLoadData('WB/DR_Gem', 'combined_model', [], true);

arMergePlotMulti(1,{'WB_DR_Dox', 'WB_DR_SN38', 'WB_DR_Gem'}, ...
   {'Dox', 'SN38', 'Gem'}, 'WB_DR');

%% qRT-PCR Data
arLoadData('qRT-PCR/gem_data', 'combined_model', [], true);
arLoadData('qRT-PCR/dox_data', 'combined_model', [], true);
arLoadData('qRT-PCR/sn38_data', 'combined_model', [], true);

arMergePlotMulti(1,{'qRT_PCR_dox_data', ...
    'qRT_PCR_sn38_data', ...
    'qRT_PCR_gem_data'}, ...
   {'Dox', 'SN38', 'Gem'}, 'qRT_PCR');

%% Compile and configure
arCompileAll;

ar.model.qPositiveX(:) = 0;
ar.model.qPlotX(:) = 0;
ar.model.qPlotU(:) = 0;

ar.config.atol = 1e-8;
ar.config.rtol = 1e-8;
ar.config.maxsteps = 1e4;
ar.config.optimizer = 1;
ar.config.optim.Display = 'off';
ar.config.optim.TolX = 1e-6;
ar.config.optim.MaxIter = 1000;
ar.config.nfine_dr_plot = 100;

arLoadPars('final_fit');
