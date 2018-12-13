%% Create *.def
arImportSBML('BIOMD0000000379', 'Tend', 100)

%% Load models & data
arInit;
ar.config.checkForNegFluxes = false;
arLoadModel('BIOMD0000000379');
arLoadData('BIOMD0000000379_data');
arCompileAll;

