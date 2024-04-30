% [pval,q] = arLRT(DeltaChi2, DeltaDF, alpha)
%
% Performs Likelihood Ratio Test
% DeltaDF   degrees of freedom
% alpha     significance level 
%
% q = TRUE: accept H0, the smaller model can not be dismissed
% q = FALSE: reject H0, dismiss smaller model in favour of larger model
% 
% Example: (ar2 ist the larger model, m1 must be nested in m2, i.e. a special case)
%       ar = ar1;
%       arCalcMerit
%       [~,m1] = arGetMerit;
%       ar = ar2;
%       arCalcMerit
%       [~,m2] = arGetMerit;
%       DeltaChi2 = m1.loglik_all - m2.loglik_all  % must be positive
%       DeltaDF = m2.npara - m1.npara; % must be positive


function [pval,q] = arLRT(DeltaChi2, DeltaDF, alpha)
if(~exist('alpha','var') || isempty(alpha))
    alpha = 0.05;
end
pval = 1-chi2cdf(DeltaChi2, DeltaDF);
q = pval > alpha;