% [Rp,V] = arGelmanRubinBrooks(d1,d2)
%
% Calculates Gelman-Rubin-Brooks diagnostic of two chains
%
%   d1       First chain, provide ar.ps samples from one runs      
%   d2       Second chain, provide ar.ps samples from different run
%            d1 and d2 must be the same length!
%
%   R       R-test value for chains
%   V       V matrix of test
%
% See also arMC3
%
% Code is minimally edited version of previously published code by 
% Benjamin Ballnus (2016): PESTO toolbox for Matlab
% Introduced in paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5482939/
% https://github.com/ICB-DCM/PESTO
%
% PESTO
%   Paul Stapor et. al (2018)
%   "PESTO: Parameter EStimation TOolbox"
%   Bioinformatics 34 (4): 705–707
%   https://doi.org/10.1093/bioinformatics/btx676 
%
% Based on Brooks Gelman 1998. Compares two multivariate data sets with
% each other using the Gelamn Rubin Brooks test.



function [Rp,V] = arGelmanRubinBrooks(d1,d2)

n1 = max(size(d1));
n2 = max(size(d2));

if n1 ~= n2
  error('Samples need the same lengths for this test.');
else
  n = n1;
end

m1 = mean(d1); 
m2 = mean(d2);
m = mean([m1;m2]);


W = 0.5/(n-1) * ((d1-repmat(m1,n,1))' * (d1-repmat(m1,n,1)) + ...
                 (d2-repmat(m2,n,1))' * (d2-repmat(m2,n,1)));
               
Bn = 0.5 * ((m1-m)'*(m1-m) + (m2-m)'*(m2-m));      

% Posterior variance-covariance estimate:
V = (n-1)/n * W + 1.5 * Bn;

% MPSRF measure (distance between matrices)
try
  lambda = max(eig(W\Bn)); % might contain small numerically caused imags and might also contain linear dependent columns
  Rp = (n-1)/n + 1.5 * lambda;  % upper bound
catch e
  Rp = inf;
  disp(e);
end

if Rp < 1     % in some cases with badly conditioned matrices, one can face negative eigenvalues.
              % in these cases we decide to assume the data sets are not
              % sufficently similar.
  Rp = inf;
end







                
               