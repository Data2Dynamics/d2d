% Upper and lower bounds are used to define a gaussian prior which has ar.mean
% at the center.
% ar.std is set as (ub-lb)/2
function arComputePriors

global ar

qunif = ar.type == 0;

means = (ar.ub + ar.lb)/2;
stds = (ar.ub - ar.lb)/2;

ar.mean(qunif) = means(qunif);
ar.std(qunif) = stds(qunif);
ar.type(qunif) = 1;
