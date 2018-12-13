% This function sets all pertubation parameters (for specifying knockout,
% knockdown and overexpression) to zero (i.e. to the wild-type setting).
% 
%   Since the bounds are also set to zero, the parameters are always in the
%   WT setting, e.g. if random initial guesses (e.g. in arFitLHS) are used.
%   This is often required/reasonable for Benchmark-Analyses.

function ResetPertubationParameters
global ar

ind = find(~cellfun(@isempty,regexp(ar.pLabel,'gene\d+_ko')));
ar.p(ind) = 0;
ar.qLog10(ind) = 0;
ar.lb(ind)=0;
ar.ub(ind)=0;
ar.qFit(ind)=0;

ind = find(~cellfun(@isempty,regexp(ar.pLabel,'sirna\d+_kd')));
ar.p(ind) = 0;
ar.qLog10(ind) = 0;
ar.lb(ind)=0;
ar.ub(ind)=0;
ar.qFit(ind)=0;

ind = find(~cellfun(@isempty,regexp(ar.pLabel,'rbs\d+_ic')));
ar.p(ind) = 0;
ar.qLog10(ind) = 0;
ar.lb(ind)=0;
ar.ub(ind)=0;
ar.qFit(ind)=0;


