% [mu, sigma] = arLogNtoN(mu_logn, sigma_logn)
%
% converts mean and sd of a log-normal distribution
% to mean and sd of the corresponding normal distribution
%	
% TODO small errors for log10 !!!


function [mu, sigma] = arLogNtoN(mu_logn, sigma_logn)

% see help lognstat
mu = log10(mu_logn^2 / sqrt(sigma_logn^2 + mu_logn^2));
sigma = sqrt(log10(sigma_logn^2 / mu_logn^2 + 1));