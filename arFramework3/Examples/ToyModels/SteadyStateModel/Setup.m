% Load models & data
arInit;
ar.config.checkForNegFluxes = false;
arLoadModel('SteadyStateModel');
arLoadData('steady',1,'csv',true);

arCompileAll;

arClearPriors('condclear');
ss_constraint_strength = 1e-3;
ar.model.condition.qSteadyState = ones(size(ar.model.condition.qSteadyState));
ar.model.condition.stdSteadyState = ss_constraint_strength*ones(size(ar.model.condition.stdSteadyState));

