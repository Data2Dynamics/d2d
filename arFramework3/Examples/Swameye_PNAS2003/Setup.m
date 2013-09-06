% Load models & data

arInit
arLoadModel('pnas_jak_stat', 1);
arLoadData('pnas_data_original', 1);
arParseModel;
arWriteCFiles;
arLink;

% load best fit parameter values
arLoadPars('BestFit');

useErrorModel = false;
if(~useErrorModel)
    % do not fit error model
    ar.config.fiterrors = -1;
    
    % show error bars instead of error model
    ar.config.ploterrors = 1;
    
    % show, but not fit input data for pEpoR
    ar.model.data.qFit(3) = 0;
    
    % fix parameter for input and error model
    arSetPars('gif_amp_sust',[],2);
    arSetPars('gif_amp_trans',[],2);
    arSetPars('gif_timescale_sust',[],2);
    arSetPars('gif_timescale_trans',[],2);
    arSetPars('scale_pEpoR',[],2);
    arSetPars('sd_pEpoR_au',[],2);
    arSetPars('sd_pSTAT_au',[],2);
    arSetPars('sd_tSTAT_au',[],2);
end

arPlot;
arPrint;
