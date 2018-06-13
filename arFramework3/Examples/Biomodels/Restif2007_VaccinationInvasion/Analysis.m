%% Create *.def
arImportSBML('BIOMD0000000294','tend',100)

%% Load models & data
Setup

%% Check
ar.config.atol=1e-10;
ar.config.rtol=1e-10;
arQplot('x')
arPlot

arCompareWithBiobaseSimulation('SIMU1448015036287.dat');

