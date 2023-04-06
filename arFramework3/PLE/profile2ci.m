% Code taken from ple.m and generalized
% 
%   x   values in horizontal direction of a profile likelihood, e.g. parameters
%       or predictions
% 
%   chi2    the profile likelihood (i.e. values in vertical direction)
%   
%   alpha   confidence level
%   df      degrees of freedom
%           [1] default: point-wise CI
% 
%   delta   threshold (alternatively to specifying alpha and df)

function CI = profile2ci(x,chi2,alpha,df,deltaChi2)
alphaOrDfGiven = ~isempty(alpha) || ~isempty(df);

if ~exist('alpha','var') || isempty(alpha)
    alpha = 1;
end

if ~exist('df','var') || isempty(df)
    df = 1;
end

if ~exist('deltaChi2','var') || isempty(deltaChi2)
    deltaChi2 = icdf('chi2',1-alpha,df);
elseif alphaOrDfGiven
    warning('profile2ci.m: Threshold and alpha/df provided: The threshold delta will overwrites alpha and df.')
end    

if length(x)~=length(chi2)
    error('length(x)~=length(chi2)')
end

% eliminate values where x is NaN
notnan = find(~isnan(x));
x = x(notnan);
chi2 = chi2(notnan);

q_chi2good = chi2<=min(chi2,[],'omitnan')+deltaChi2 & ~isnan(chi2);

if(min(x(q_chi2good))==min(x(~isnan(chi2))))
    CI = -Inf;
    % bounds are infinite if chi2 does not exceed boundary
else
    kind = find(x==min(x(q_chi2good)));
    try
        CI = interp1(chi2([kind kind-1]), x([kind kind-1]), min(chi2)+deltaChi2);
        % interpolate to increase accuracy of confidence interval
    catch
        CI = NaN;  % e.g. isinf(chi2(kind-1))
    end
end
if(max(x(q_chi2good))==max(x(~isnan(chi2))))
    CI(2) = Inf;
else
    kind = find(x==max(x(q_chi2good)));
    try
        CI(2) = interp1(chi2([kind kind+1]), x([kind kind+1]), min(chi2)+deltaChi2);
    catch
        CI(2) = NaN;
    end

end
