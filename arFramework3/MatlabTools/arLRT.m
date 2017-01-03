% LRT
%
% [pval,q] = arLRT(DeltaChi2, DeltaDF, alpha)
% pval = 1-chi2cdf(DeltaChi2, DeltaDF);
% q = pval > alpha;
%
% q = TRUE: accept H0, the smaller model can not be dismissed
% q = FALSE: reject H0, dismiss smaller model in favour of larger model

function [pval,q] = arLRT(DeltaChi2, DeltaDF, alpha)
pval = 1-chi2cdf(DeltaChi2, DeltaDF);
q = pval > alpha;