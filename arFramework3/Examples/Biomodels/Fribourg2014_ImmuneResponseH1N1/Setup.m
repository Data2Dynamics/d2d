%% Create *.def
arImportSBML('BIOMD0000000528_d2d','tend',100)

%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000528_d2d');
arLoadData('BIOMD0000000528_d2d_data');
arCompileAll;

ar.config.atol=1e-8;
ar.config.rtol=1e-8;

