%% 
% This file describes the use of the  multistart fit merging approach 
% based on the Nudged Elastic Band (NEB) path finding method using the
% Swameye et al. (2003) model with 10 fitted paramters and enables to 
% generate the results from Tönsing et al. (2019)


%% Multistart Fit of Swameye model
close all, clear all

arInit
arLoadModel('pnas_jak_stat');
arLoadData('pnas_data_original');
arCompileAll;

% load best fit parameter values
arLoadPars('BestFit');

% do not fit input splines
ar.qFit(13:17) = 0; % sp's input splines

% do not fit pEpo error
ar.qFit(10) = 0; % sd_pEpoR_au

ar.config.useFitErrorCorrection = 0; % Turn Bessel Correction off
ar.config.fiterrors = 1; % Use error model of the Swameye model

ar.config.optim.TolX = 1e-1; % intentionally suboptimal optimizer setting
%ar.config.optim.TolX = 1e-6; % standard optimizer setting

arFit

%Perform multistart fitting sequence with 50 fits and with random seed '1336'
arFitLHS(50,1336)

% Save Results for later use in NEB approach
arSave('Multistart_Results')


%% Merging Fits using NEB method

% *** PREREQUISITES ***
%
% 1) prepare NEB model file template: e.g. Models/pnas_jak_stat_NEB_XXX.def
%
%
% All NEB parameters have to be defined in the first lines of the template
% You may select the fitted parameters by the ar.qFit in the model a 
% generate a list in the appropriate format by of using the function

arNEBdefineList 

% afterwards. The output needs to be copied to the beginning of the 
% NEB model file template (model file with suffix '_NEB_XXX') 
% 
% (if applicable, the same procedure has to be done for data def files)
%
%
% 2) generate explicit NEB node model files from template(s)
%
% If the template, e.g. Models/pnas_jak_stat_NEB_XXX.def exists,
% the function 

arNEBMakeModelFiles('pnas_jak_stat',25) 

% generates the 25 model files for the 25 free NEB nodes
%
% (if applicable, the same procedure has to be done for data def files,
% using arNEBMakeDataFilesCSV, arNEBMakeDataFilesXLS or arNEBMakeDataFilesXLSX)



% Initialize NEB band using multiple model instances

close all, clear all

arInit
arLoadModel('pnas_jak_stat');
arLoadData('pnas_data_original',1);

for inode = 1:25
    arLoadModel(['pnas_jak_stat_NEB_' sprintf('%03d', inode) ]);
    arLoadData('pnas_data_original',inode+1);
 
end
arLoadModel('pnas_jak_stat_NEB_end');
    arLoadData('pnas_data_original',inode+2);

arCompileAll;

% Save compiled NEB band for later use
arSave('NEB_25steps')


 
% Load multistart results (waterfall plot) 
arLoadOnlyMultistartResults % select respective multistart results from above

ar.config.useFitErrorCorrection = 0; % Turn Bessel Correction off
ar.config.fiterrors = 1;% Use error model of the Swameye model

ar.config.optim.TolX = 1e-1; % intentionally suboptimal optimizer setting
%ar.config.optim.TolX = 1e-6; % standard optimizer setting



% Initialize NEB path finding and merging
arInitNEB

% choose spring constants to be tested
springs = [200,100,50,20,10,5,2,1];

qsaveplots = 0; % turn saving of plots off

% Calculate NEBs and merge fits #1-#50 of waterfall plot
[merged_lhs] = arNEBMergeLHS(1:50, springs, qsaveplots);

% Plot results (merged waterfall plot)
arNEBPlotMerged 

% Save results
arSave('NEB_Results')

