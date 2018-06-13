%% source: DDmore http://repository.ddmore.foundation/model/DDMODEL00000103#Overview
% data103_IDall: data file from DDmore, 569 individuals, a lot with <5 mesurements,
%   takes very long to compute, only one individual parameter possible
% data103_IDlong: data of all individuals with at least 5 mesurements
% data103_IDshort: data of idividual 508-523, similar results to the long file,
%   significant shorter computation time

% initialize model
arInit
arLoadModel('mod103');
arLoadData('data103_IDshort');
arCompileAll;

% changing the prior for the standard deviation of the individual parameters
ind = find(~cellfun(@isempty,regexp(ar.pLabel,'sd')));
ar.type(ind) = 1;
ar.mean(ind) = 0.5;
ar.std(ind) = 0.3;
ar.lb(ind) = -0.3;
ar.ub(ind) = 1;
ar.p(ind) = 0;
ar.p(6) = -5;
ar.qFit(6) = 0;


