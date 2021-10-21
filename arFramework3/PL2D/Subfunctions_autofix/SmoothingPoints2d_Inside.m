function [indxs,indys,dirs,margin] = SmoothingPoints2d_Inside(output_msg)

% Proposes 2d-profile discontinuities. This is done by finding the largest
% chi2-difference between the common grid points of two neighboring profiles 
% and checking whether this is unusually large by comparing the difference
% to the difference to the previous and next prediction value at the same
% parameter index indx.
%
% indxs:   Parameter index in ar.ple2d.smooth of discontinuous point
% indys:   Prediction index in ar.ple2d.smooth of discontinuous point
% dirs:    Direction in which correction of profile is necessary
% margin:  Chi2-difference between profiles of possible discontinuity
%
% This function is needed if the behavior of the profiles is continuous at
% the bounds, but not between the bounds.
%
% See also: rmProfileJump2d , SmoothingPoints2d_Excess,
% SmoothingPoints2d_Bounds

if ~exist('output_msg','var') || isempty(output_msg)
    output_msg = 1;
end

global ar

n = length(ar.ple2d.raw.predsteps);
jumpfactor = ar.ple2d.config.autofix.jumpfactor;

smooth2d;

indxs  = NaN(1,n-1);
indys  = NaN(1,n-1);
dirs   = NaN(1,n-1);
margin = NaN(1,n-1);

for indy_tmp = 1:(n-1)
    
    [margin_tmp,indxmax_tmp] = max(abs(ar.ple2d.smooth.zq(indy_tmp+1,:)...
        - ar.ple2d.smooth.zq(indy_tmp,:)));
    sig = sign(ar.ple2d.smooth.zq(indy_tmp+1,indxmax_tmp) ...
        - ar.ple2d.smooth.zq(indy_tmp,indxmax_tmp));
    % sig contains information about direction
    
    % Check whether maximal difference between profiles is sufficiently
    % large. If not, iteration skips this value.
    if (margin_tmp < ar.ple2d.config.autofix.eps_inside) || ...
            isnan(margin_tmp)
        continue
    else
        % If difference is larger (by a certain factor) than previous
        % and next difference, this might be a jump:
        if indy_tmp == 1 
            % Edge values are exception
            diffs_tmp = abs(diff(ar.ple2d.smooth.zq(1:3,indxmax_tmp)));
            if jumpfactor*diffs_tmp(2) > diffs_tmp(1)
                continue
            end
        elseif indy_tmp == n-1
            diffs_tmp = abs(diff(ar.ple2d.smooth.zq((n-2):n,indxmax_tmp)));
            if jumpfactor*diffs_tmp(1) > diffs_tmp(2)
                continue
            end
        else
            diffs_tmp = abs(diff(ar.ple2d.smooth.zq((indy_tmp-1):(indy_tmp+2),...
                indxmax_tmp)));
            if (jumpfactor*diffs_tmp(1) > diffs_tmp(2)) || ...
                    (jumpfactor*diffs_tmp(3) > diffs_tmp(2)) 
                continue
            end
        end
    end
    
    % If iteration reaches this point, a jump was detected. Save this
    % value:
    margin(indy_tmp) = margin_tmp;
    indxs(indy_tmp) = indxmax_tmp;
    if sig > 0
        indys(indy_tmp) = indy_tmp;
        dirs(indy_tmp) = 1;
    elseif sig < 0
        indys(indy_tmp) = indy_tmp+1;
        dirs(indy_tmp) = -1;
    end
    
end

q_tmp = ~isnan(indxs);
indxs = indxs(q_tmp);
indys = indys(q_tmp);
dirs  = dirs(q_tmp);
margin = margin(q_tmp);


if output_msg == 1
    if ~isempty(indxs)
        fprintf(['\n %i potential profile smoothing points have been',...
            ' found (inside check). \n'],length(indxs));
    else
        fprintf(['\n No potential profile smoothing points have been',...
            ' found (inside check). \n']);
    end
end

end