function arComputePriors

global ar

qunif = ar.type == 0;

means = (ar.ub + ar.lb)/2;
stds = (ar.ub - ar.lb)/2;

ar.mean(qunif) = means(qunif);
ar.std(qunif) = stds(qunif);
ar.type(qunif) = 1;
