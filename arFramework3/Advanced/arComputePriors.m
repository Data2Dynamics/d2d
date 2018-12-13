% arComputePriors
% 
% Upper and lower bounds are used to define a gaussian priors.

% ar.mean is set as the average of upper bounds ar.ub and lower bounds ar.lb
% ar.std is set as (ar.ub-ar.lb)/2
% ar.type is set to 1 which indicates a Gaussian prior

function arComputePriors

global ar

qunif = ar.type == 0;

means = (ar.ub + ar.lb)/2;
stds = (ar.ub - ar.lb)/2;

ar.mean(qunif) = means(qunif);
ar.std(qunif) = stds(qunif);
ar.type(qunif) = 1;
