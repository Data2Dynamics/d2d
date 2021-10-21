function [indxs,indys,dirs,margin] = SmoothingPoints2d_Excess(lbub,output_msg)

% [indys,dirs,indxs,margin] = SmoothingPoints2d_Excess(lbub,output_msg)
%
% This function lists all lower (lbub = -1) or upper (lbub = 1) parameter 
% bound values and finds possible jumps depending on unusual 
% (i.e. discontinuous) behavior of the bounds. If this is the case, one 
% profile extends unusually much farther than the other, thus, this is the 
% excess method. 
%
% Computationally cheap.
%
% Returns the prediction indices (indys) of the starting profiles for the 
% extension and the prediction directions (dirs) in which parameter profiles 
% are extended. The smooth grid parameter indices indxs give the parameter from
% which extension should be tried, chosen by having the minimal chi^2 value
% in the profile excess.
%
% indxs:   Parameter index in ar.ple2d.smooth of discontinuous point
% indys:   Prediction index in ar.ple2d.smooth of discontinuous point
% dirs:    Direction in which correction of profile is necessary
% margin:  Chi2-margin by which the discontinuity was detected
%
% See function comments for additional information.
%
% See also: rmProfileJumps2d , SmoothingPoints2d_Bounds,
% SmoothingPoints2d_Inside

if ~exist('lbub','var') || isempty(lbub)
    fprintf(['\n ERROR SmoothingPoints2d_Excess: \n Please specify',...
        ' whether lower (lbub = -1) or upper (lbub = 1) is of interest. \n']);
    return
end
if ~exist('output_msg','var') || isempty(output_msg)
    output_msg = 1;
end

global ar

jumpfactor = ar.ple2d.config.autofix.jumpfactor;
scalefactor = ar.ple2d.config.autofix.scalefactor;

smooth2d;

n = length(ar.ple2d.raw.predsteps);
thresh = scalefactor*(max(ar.ple2d.raw.plpar,[],'all') ...
    - min(ar.ple2d.raw.plpar,[],'all'));

% Lower or Upper bound:
if lbub == -1
    bounds = min(ar.ple2d.raw.plpar);
    txt = 'lower';
elseif lbub == 1
    bounds = max(ar.ple2d.raw.plpar);
    txt = 'upper';
end

% Find jump postions (prediction index) but not direction.
%
% Criterion: For each difference of bound values compare the two neighboring
% differences. If it is larger than both of its neighbors by a pre-defined
% factor, this is likely a discontinuity.
indys = [1,2,n-1,n];
dirs = [1,-1,1,-1];
diff_bounds = abs(diff(bounds)); % absolute value!
for ii = 2:(n-2)
    if (jumpfactor*diff_bounds(ii-1) < diff_bounds(ii)) && ...
            (jumpfactor*diff_bounds(ii+1) < diff_bounds(ii))
        indys = [indys,ii,ii+1];
        dirs = [dirs,1,-1];
        % Only knowning the difference leaves ambiguity
    end
end

% Sort out the right index and direction:
%
% Criterion: Find which of profiles extends further.
q_tmp = logical(zeros(size(indys)));
for ii = 1:length(indys)
    diff_tmp = lbub*(bounds(indys(ii))-bounds(indys(ii)+dirs(ii)));
    if (diff_tmp > 0) && (diff_tmp > thresh)
        % thresh gives an absolute limit on how large the difference must
        % be at least.
        q_tmp(ii) = 1;
    end
end
indys = indys(q_tmp);
dirs = dirs(q_tmp);

% Determine final trial points.
%
% Criterion: For the starting profile and the profile next to it (in the
% right direction), by construction one profile extends further to one
% side than the other. Check if any of these values is under the confidence
% threshold (extendend by ar.ple2d.config.autofix.excess_factor).
% If this is the case, the neighboring profile will likely also have
% values under the threshold there which were not found originally.
q_grid = ~isnan(ar.ple2d.smooth.zq);
indxs = NaN(size(indys));
margin = NaN(size(indys));
for ii = 1:length(indys)
    indy_now = indys(ii);
    indy_nxt = indys(ii) + dirs(ii);
    if lbub == -1
        indx_now = find(q_grid(indy_now,:),1,'first');
        indx_nxt = find(q_grid(indy_nxt,:),1,'first');
    elseif lbub == 1
        indx_now = find(q_grid(indy_now,:),1,'last');
        indx_nxt = find(q_grid(indy_nxt,:),1,'last');
    end
    
    if indx_now - indx_nxt == 0
        continue
    else
        indx_list = (indx_nxt+lbub):lbub:indx_now;
        % For these parameter indices there is a profile excess
        chi2_tmp = ar.ple2d.smooth.zq(indy_now,:);
        [min_excess,subindx_tmp] = min(chi2_tmp(indx_list));
        indx_tmp = indx_list(subindx_tmp);
        % Lowest chi of all excess values
        margin(ii) = ar.ple2d.config.autofix.excess_factor *...
            icdf('chi2',ar.ple2d.config.bounds.alpha_bounds,1)...
            - (min_excess - min(chi2_tmp));
        if margin(ii) <= 0
            continue
        else
            indxs(ii) = indx_tmp;
        end
        % Under threshold? If yes, save value.
    end
end
q_tmp = ~isnan(indxs);
indxs = indxs(q_tmp);
indys = indys(q_tmp);
dirs  = dirs(q_tmp);
margin = margin(q_tmp);

if output_msg == 1
    if ~isempty(indxs)
        fprintf(['\n %i potential %s bound profile smoothing points have been',...
            ' found (excess check). \n'],length(indxs),txt);
    else
        fprintf(['\n No potential %s bound profile smoothing points have been',...
            ' found (excess check). \n'],txt);
    end
end

end

