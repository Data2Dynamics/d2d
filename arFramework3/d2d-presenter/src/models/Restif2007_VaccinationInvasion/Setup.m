%% Create *.def
arImportSBML('BIOMD0000000294',100)

%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000294');
arLoadData('BIOMD0000000294_data');
arCompileAll;

%% Check
ar.config.atol=1e-10;
ar.config.rtol=1e-10;
%arQplot('x')
%arPlot

%arCompareWithBiobaseSimulation('SIMU1448015036287.dat');

