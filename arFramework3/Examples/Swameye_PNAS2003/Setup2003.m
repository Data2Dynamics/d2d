% Old version of the Swameye 2003 model with linear interpolation

% Load models & data
arInit
arLoadModel('pnas_jak_stat3');
arLoadData('pnas_data_original');
arCompileAll;

% load best fit parameter values
arLoadPars('BestFit');

% do not fit error model
ar.config.fiterrors = -1;
    
% show error bars instead of error model
ar.config.ploterrors = 1;
    
% show, but not fit input data for pEpoR
ar.model.data.qFit(3) = 0;
    
% fix parameter for input and error model
arSetPars('scale_pEpoR',[],2);
arSetPars('sd_pEpoR_au',[],2);
arSetPars('sd_pSTAT_au',[],2);
arSetPars('sd_tSTAT_au',[],2);

arFit;

arPlot;
arPrint;
