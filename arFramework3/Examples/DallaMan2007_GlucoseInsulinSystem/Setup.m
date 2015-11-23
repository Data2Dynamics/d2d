%% Create *.def
arImportSBML('BIOMD0000000379',100)

%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000379');
arLoadData('BIOMD0000000379_data');
arCompileAll;

%% Check
ar.config.atol=1e-2;  % Attention: 1e-10 leads to wrong integration
ar.config.rtol=1e-2;  % Attention: 1e-10 leads to wrong integration
arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1448030646732.dat');

