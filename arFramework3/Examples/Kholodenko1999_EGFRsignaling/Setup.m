%%
% addpath('E:\clemens\systemBiologie_cvs\Programmieren\Matlab\matlab-pathway\Library\libSBML-5.12.0-matlab')
% arImportSBML('BIOMD0000000048',100)
arImportSBML('BIOMD0000000048_d2d',100)

%%
arInit;
ar.config.checkForNegFluxes = false
% arLoadModel('BIOMD0000000048');
% arLoadData('BIOMD0000000048_data');
arLoadModel('BIOMD0000000048_d2d');
arLoadData('BIOMD0000000048_d2d_data');
arCompileAll;


%% Comparison with Biobase:
ar.config.atol=1e-14;
ar.config.rtol=1e-14;

arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1447940350922.dat');
