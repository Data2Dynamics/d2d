%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000009_d2d');
arLoadData('BIOMD0000000009_d2d_data');
arCompileAll;

ar.config.atol=1e-10;
ar.config.rtol=1e-10;
