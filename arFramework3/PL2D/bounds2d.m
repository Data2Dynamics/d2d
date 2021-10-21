function [pred,alpha_pred,bounds,par_min] = bounds2d(alphapar_max)
% [pred,alpha_pred,bounds,par_min] = bounds2d(alphapar_max)
%
% Returns the parameter bounds of profiles with different valdiation
% confidence levels which are controlled in ar.ple2d.config.bounds. The
% corresponding parameter profiles are constructed and saved by the
% function findvalidationpoints.
%
% Input:
% alphapar_max = Parameter profile confidence level for which confidence
%                   bounds are calculated [ar.ple2d.config.bounds.alpha_bounds] 
%
% Output:
% pred: Sampled prediction points 
% alpha_pred: Validation confidence levels corresponding to predictions in pred
% bounds: Matrix of boundaries for all parameter profiles corresponding to
%           each prediction in pred
% par_min: Vector of minimal parameter values for each parameter profile.
%            Caution: This may very well be another parameter value than
%            for the independent validation profile calculation, if the
%            validation profile has not found the minimum for this prediction.
%
% Use fields ar.ple2d.config.bounds.vpl_mode to control whether original or grid
% parameter profiles should be used and ar.ple2d.config.bounds.save_mode to
% control whether output should be saved for plotting and analysis.
%
% See also: gen2d, score2d

global ar

% This argument is necessary because it is repeatedly changed in score2d 
if ~exist('alphapar_max','var') || isempty(alphapar_max)
    alphapar_max = ar.ple2d.config.bounds.alpha_bounds;
end

%% Find averaging points and corresponding profiles

% For the score calculation findvalidationpoints is calculated
% independently and save_mode is set to 1, such that this
% findvalidationpoints is only called once. When explicitly calling bounds2d,
% findvalidationpoints is always reevaluated since inputs may have changed.
if ar.ple2d.config.bounds.save_mode == 2
    findvalidationpoints;
end
pred = ar.ple2d.bounds.pred_points;
alpha_pred = ar.ple2d.bounds.alpha_pred;
par = ar.ple2d.bounds.pars_profile;
chi2 = ar.ple2d.bounds.chi2s_profile;

%% Find boundaries of each profile
% Idea: Loop over every prediction (for which the respective parameter profile 
% was saved) and find the intersections of the profile with a specified
% confidence level

eps = 10^-1;
bounds = NaN(length(pred),4); %is extended if necessary
q_lb = logical(zeros(length(pred),1));
q_ub = logical(zeros(length(pred),1));
par_min = NaN(length(pred),1);
for ii = 1:length(pred)
    chi2_profile = chi2(:,ii)';
    [min_chi2,indmin_chi2] = min(chi2_profile);
    % bounds are always calculated wrt to the actual minimum, not to the one 
    % given from the independent validation profile.
    chi2_profile = chi2_profile(~isnan(chi2_profile));
    chi2_profile = chi2_profile-min_chi2;   
    par_profile = par(:,ii)';
    par_min(ii) = par_profile(indmin_chi2);
    par_profile = par_profile(~isnan(par_profile));
    intersects = findsamplepoints(par_profile,chi2_profile,...
        icdf('chi2',alphapar_max,1));
    % Vector of intersections of parameter profile and confidence threshold
    
    if length(intersects) > size(bounds,2)
        bounds = [bounds,NaN(size(bounds,1),length(intersects) - size(bounds,2))];
    end
    bounds(ii,1:length(intersects)) = intersects;
    % Insert values into bounds (which is extended if it is not big enough)

    % Use hard parameter bounds if profile does not reach confidence level,
    % append them to bounds later
    if (abs(min(par_profile)-ar.lb(ar.ple2d.general.idpar)) < eps)  && ...
            (chi2_profile(1) < icdf('chi2',alphapar_max,1))
        q_lb(ii) = true;
    end
    if (abs(max(par_profile)-ar.ub(ar.ple2d.general.idpar)) < eps) && ...
            (chi2_profile(end) < icdf('chi2',alphapar_max,1))
        q_ub(ii) = true;
    end

    % There are no reasonable bounds for an uneven amount of
    % intersections if no hard boundary is reached:        
    if (rem(length(intersects),2) == 1) && ~q_ub(ii) && ~q_lb(ii)
        fprintf(['\n ERROR bounds2d: A Profile intersects the confidence',...
            ' level an uneven number of times \n without reaching the',...
            ' hard parameter bounds. Try using autofix2dprofile to extend',...
            ' profiles or try setting ar.ple2d.config.bounds.vpl_mode',...
            ' = 1.\n']);
        return
    end

end

% Append hard boundaries to bounds:
bounds = [NaN(size(bounds,1),1),bounds,NaN(size(bounds,1),1)];
bounds(q_lb,1) = ar.lb(ar.ple2d.general.idpar);
bounds(q_ub,end) = ar.ub(ar.ple2d.general.idpar);

if ar.ple2d.config.bounds.save_mode == 2
    ar.ple2d.bounds.bounds = bounds;
    ar.ple2d.bounds.alpha_bounds = alphapar_max;
end

end

