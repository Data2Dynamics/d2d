arImportSBML('BIOMD0000000012',1000)

%%
arInit;
arLoadModel('BIOMD0000000012');
arLoadData('BIOMD0000000012_data');
arCompileAll;


%%
ar.config.atol=1e-10;
ar.config.rtol=1e-10;
arQplot('x')
arPlot

%%
arCompareWithBiobaseSimulation('SIMU1447941514719.dat');