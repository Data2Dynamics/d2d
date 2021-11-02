function [indxs,indys,dirs,margin] = SmoothingPoints2d_Bounds(lbub,output_msg)

% [indxs,indys,dirs,margin] = SmoothingPoints2d_Bounds(lbub,output_msg)
%
% Finds discontinuities at the boundaries of the 2D-profile by proposing
% parameter steps not in parameter direction, but in data direction. If
% such a proposed step leads to an improved 2dprofile value, there is a
% discontinuity. 
%
% lbub: 1 for upper bound, -1 for lower bound 
% 
% indxs:   Parameter index in ar.ple2d.smooth of discontinuous point
% indys:   Prediction index in ar.ple2d.smooth of discontinuous point
% dirs:    Direction in which correction of profile is necessary
% margin:  Improvement of 2dprofile value over existing value
%
% This function is in its current state not particularly efficient with
% respect to computational cost, since the amount of local fits performed 
% is about the number of sampled profiles. However, all discontinuities
% found are actually real and not just possibly discontinuous.
% 
% See also: rmProfileJump2d , SmoothingPoints2d_Excess,
% SmoothingPoints2d_Inside

if ~exist('lbub','var') || isempty(lbub)
    fprintf(['\n ERROR SmoothingPoints2d_Bounds: \n Please specify',...
        ' whether lower (lbub = -1) or upper (lbub = 1) is of interest. \n']);
    return
end
if ~exist('output_msg','var') || isempty(output_msg)
    output_msg = 1;
end

global ar

ar_old = arDeepCopy(ar);

try
    
    smooth2d; %Generate a grid with regular data
    q_zq = ~isnan(ar.ple2d.smooth.zq);
    q_trial = zeros(size(q_zq));
    % q_trial contains information about direction. q_trial = 2 means that
    % positive and negative direction should be checked.
    
    for ii = 2:(size(q_zq,1)-1)
        % Loop over every parameter profile and find parameter boundaries
        if lbub == -1
            mid = find(q_zq(ii,:) == 1,1,'first');
            if mid == 1
                continue
            end
            txt = 'lower';
        elseif lbub == 1
            mid = find(q_zq(ii,:) == 1,1,'last');
            if mid == size(q_zq,2)
                continue
            end
            txt = 'upper';
        end
        
        % Idea: If there is a sampled point on the grid above or below (in
        % data direction) the boundary point, make a step from there to
        % this boundary point and check whether chi^2 changes.        
        if (q_zq(ii-1,mid) == 1) && (q_zq(ii+1,mid) == 1)
            q_trial(ii,mid) = 2;
        elseif q_zq(ii-1,mid) == 1
            q_trial(ii,mid) = -1;
        elseif q_zq(ii+1,mid) == 1
            q_trial(ii,mid) = 1;
        end
    end
    
    [m,d,idpred,iddata] = Set2dDataPlaceholder;
    
    inds_trial = find(q_trial); %all points which must be checked
    dirs = NaN(2*length(inds_trial),1);
    indxs = NaN(2*length(inds_trial),1);
    indys = NaN(2*length(inds_trial),1);
    margin = NaN(2*length(inds_trial),1);
    
    fprintf(['\n Look for potential bound smoothing points',...
        ' in %i instances ... \n'],length(inds_trial))
    
    for ii = 1:length(inds_trial)
        
        % Indices of point which is currently checked
        [indy_tmp,indx_tmp] = ind2sub(size(q_zq),inds_trial(ii));
        
        % Allow for two directions to be checked
        if q_trial(inds_trial(ii)) ~= 2
            direc = q_trial(inds_trial(ii));
        else
            direc = [-1,1];
        end
        
        for jj = 1:length(direc)
            
            % Index of the point whichs pecifies the new trial parameter
            indy_prev = indy_tmp + direc(jj);
            
            % Take trial parameter from the previuosly calculated
            % parameter profile:
            xs = ar.ple2d.raw.plpar(:,indy_prev);
            ps = ar.ple2d.raw.par{indy_prev};
            ps = ps(~isnan(xs),:);
            xs = xs(~isnan(xs));
            p_trial = interp1(xs,ps,ar.ple2d.smooth.xq(1,indx_tmp));
            ar.p = p_trial;
            
            % Set data point and make the fit. Note the modifications necessary
            % to get the right objective functon value.
            ar.model(m).data(d).yExp(iddata,idpred) = ...
                ar.ple2d.smooth.yq(indy_tmp,1);
            chi2 = UsePLEFit(ar.ple2d.general.idpar);
            
            % Developer's note:
            % ar.ple2d.smooth.zq(indy_tmp,indx_tmp) + ar.ple2d.raw.chi2_min + ...
            % ar.ple2d.smooth.offset-chi2s(ii) should be more or less exactly
            % zero, but it is not?
            
            % Check whether this trial parameter finds a better optimum:
            chi2 = chi2 - ar.ple2d.raw.chi2_min;
            margin_tmp = ar.ple2d.smooth.zq(indy_tmp,indx_tmp) - chi2;
            if margin_tmp > ar.ple2d.config.autofix.eps_outside
                indxs(ii-1+jj) = indx_tmp;
                indys(ii-1+jj) = indy_prev;
                dirs(ii-1+jj) = -direc(jj);
                margin(ii-1+jj) = margin_tmp;
            end
            
        end
    end
    
catch exception
    fprintf(['ERROR SmoothingPoints2d_Bounds: Resetting ar struct.',...
        ' Error message: \n %s \n Line: %s \n'] ,...
        exception.message, sprintf('%i, ',exception.stack.line));
    ar = ar_old;
    return
end

ar = ar_old;
indxs = indxs(~isnan(indxs));
indys = indys(~isnan(indys));
dirs = dirs(~isnan(dirs));
margin = margin(~isnan(margin));

if output_msg == 1
    if ~isempty(indxs)
        fprintf(['\n %i %s bound profile smoothing points have been confirmed.',...
            ' (bounds check) \n'],length(indxs),txt);
    else
        fprintf(['\n No %s bound profile smoothing points have been found.',...
            ' (bounds check) \n'],txt);
    end
end
end

