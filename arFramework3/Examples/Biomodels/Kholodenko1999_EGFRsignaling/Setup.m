%%
arInit;
ar.config.checkForNegFluxes = false
% arLoadModel('BIOMD0000000048');
% arLoadData('BIOMD0000000048_data');
arLoadModel('BIOMD0000000048_d2d');
arLoadData('BIOMD0000000048_d2d_data');
arCompileAll;

ar.config.atol=1e-14;
ar.config.rtol=1e-14;

