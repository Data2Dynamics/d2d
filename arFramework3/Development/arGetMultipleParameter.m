% function [ parameter ] = arGetMultipleParameter([n], [pl])
%
% After multi-start fitting this function chooses the best resulting fits.
% 
%   pl          logical, that indicates whether plots should be produced [false]
%   n           1) integer n>1 chooses the n best fits according to ar.chi2s
%               2) 0<n<1 interprets n as a quantile and chooses those fits that
%                  are within the thresholds
% 
%   parameter   the parameter sets that are chosen by the above criterion.
% 
% If pl is set to true, arPlotMulti is called and plots the observables of
% these parameter sets.
% 
% See also arPlotMulti arFitLHS arFits

function [ parameter ] = arGetMultipleParameter(n, pl)
%
%
% pl   [false]
% n    10

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('n','var'))
	n = 10;
end

if(~exist('pl','var'))
	pl = false;
end

if(~exist('saveToFile','var'))
	saveToFile = false;
end

if(n < 1) %interpret n as quantil
    quant = chi2inv(n, numel(ar.p));
    maxChi2 = arGetMerit('chi2fit') + quant;
    indices = find(ar.chi2s < maxChi2);
    parameter = ar.ps(indices,:);  
    paraSize = size(parameter);
    fprintf('%i parameter sets selected \n',paraSize(1));
else % n as number
    if(numel(ar.chi2s)-sum(isnan(ar.chi2s)) < n)
        fprintf('number of requested parameter sets is greater as the number of available parameter sets \n')
        n = numel(ar.chi2s)-sum(isnan(ar.chi2s));
        fprintf('setting number of parameter sets to %i \n', n)
    end    
        
    sorted = sort(ar.chi2s);
    indices = zeros([1 n]);
    for i = 1:n
        indices(i) = find(ar.chi2s == sorted(i));
    end

    parameter = ar.ps(indices,:);
end

if(pl)
    arPlotMulti(ar.ps(indices,:)); 
end


end

