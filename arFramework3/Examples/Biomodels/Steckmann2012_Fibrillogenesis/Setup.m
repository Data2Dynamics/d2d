%%
% addpath('E:\clemens\systemBiologie_cvs\Programmieren\Matlab\matlab-pathway\Library\libSBML-5.12.0-matlab')
arImportSBML('BIOMD0000000533','tend',100)

%%
arInit;
ar.config.checkForNegFluxes = false;
arLoadModel('BIOMD0000000533');
arLoadData('BIOMD0000000533_data');
arCompileAll;


%% Comparison with Biobase:
ar.config.atol=1e-10;
ar.config.rtol=1e-10;

arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1447947782219.dat');
