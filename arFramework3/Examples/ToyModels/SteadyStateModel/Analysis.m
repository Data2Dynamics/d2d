% Load models & data
Setup


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