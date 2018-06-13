%% Load models & data
arInit
ar.config.checkForNegFluxes = false
arLoadModel('BIOMD0000000294');
arLoadData('BIOMD0000000294_data');
arCompileAll;

ar.config.atol=1e-10;
ar.config.rtol=1e-10;
