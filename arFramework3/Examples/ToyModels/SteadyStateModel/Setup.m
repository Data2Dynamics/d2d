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

global ar;

tic;
for c = 1 : 5
ar.p=15*ones(size(ar.p));
% Simulate the model dynamically
ar.config.skipSim = 0;
arFit;
ar.config.skipSim = 0;
end
dyn=toc;

tic;
for c = 1 : 5
% Don't actually simulate the model dynamically
ar.p=15*ones(size(ar.p));
ar.config.skipSim = 1;
arFit;
ar.config.skipSim = 0;
end
nondyn=toc;

fprintf( 'Dynamic simulation %g seconds, only steady state %g seconds\n', dyn, nondyn );