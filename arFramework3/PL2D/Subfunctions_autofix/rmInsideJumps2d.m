function [count_improvements,count_trial] = rmInsideJumps2d

% [count_improvements,count_trial] = rmInsideJumps2d
%
% Removes discontinuities in the 2d profile due to non-converged fits. They
% are found by checking whether they exceed both their neighboring points with
% same profile parameter but other prediction. It is then attempted to fit
% with the starting parameter vector set to the one of the neighboring
% point. If the 2dprofile value improved, take the better value.

global ar

ar_old = arDeepCopy(ar);

try
    % Temporaily store results in dummy struct:
    dummy.chi2 = ar.ple2d.raw.chi2;
    dummy.plpar = ar.ple2d.raw.plpar;
    dummy.par = ar.ple2d.raw.par;
    
    % Generate regular 2D-profile grid:
    smooth2d;
    
    % Find trial points for smoothing. If too many points are found,
    % increase detection limit:
    q_trial = NaN(10^2+1,1);
    eps = ar.ple2d.config.autofix.eps_sample/2;
    while size(q_trial,1) > 10^2
        q_trial = FindInsideJumps(eps);
        eps = 2*eps;
    end
    
    if ~isempty(q_trial)
        fprintf(['\n Found %i potential discontinuities with threshold %0.4g.',...
            ' Attempting to remove them... \n'],size(q_trial,1),eps);
    else
        fprintf('\n No discontinuities have been detected. \n')
    end
    
    % Initialize data point:
    [m,d,idpred,iddata] = Set2dDataPlaceholder;
    
    % If smoothing candidates have been found, begin smoothing:
    count_improvements = 0;
    count_trial = size(q_trial,1);
    
    if ~isempty(q_trial)
        for ii = 1:size(q_trial,1)
            for jj = 1:2
                % Attempt smoothing from each data-direction:
                if jj == 1
                    indy_trial = q_trial(ii,2);
                    direc = 1;
                else
                    indy_trial = q_trial(ii,3);
                    direc = -1;
                end
                
                % Set right data point value:
                ar.model(m).data(d).yExp(iddata,idpred) = ...
                    ar.ple2d.smooth.yq(indy_trial+direc,1);
                
                % Initialize parameter vector by virtue of parameter
                % profile above/below:
                ptrial_tmp = ptrialUpdate(dummy,indy_trial,indy_trial +direc,...
                    q_trial(ii,1),[]);
                if isempty(ptrial_tmp)
                    continue
                end
                
                % Optimize with new parameter vector:
                ar.p = ptrial_tmp;
                chi2 = UsePLEFit(ar.ple2d.general.idpar);
                chi2_new = chi2 - ar.ple2d.raw.chi2_min;
                
                [~,indx_trial] = min(abs(...
                    dummy.plpar(:,indy_trial+direc) - q_trial(ii,1)));
                % This is just to find the vector of raw profile parameter
                % values. The actual difference should be zero by construction.
                
                % If profile value improved, save the new values:
                if chi2_new + 0.01 < dummy.chi2(indx_trial,indy_trial+direc)
                    dummy.chi2(indx_trial,indy_trial+direc) = chi2_new;
                    dummy.plpar(indx_trial,indy_trial+direc) = ...
                        ar.p(ar.ple2d.general.idpar);
                    dummy.par{indy_trial+direc}(indx_trial,:) = ar.p;
                    count_improvements = count_improvements +1;
                end
            end
        end
    end
    
catch exception
    fprintf(['ERROR rmInsideJumps2d: Resetting ar struct.',...
        ' Error message: \n %s \n Line: %s \n'] ,...
        exception.message, sprintf('%i, ',exception.stack.line));
end

ar = ar_old;

ar.ple2d.raw.chi2 = dummy.chi2;
ar.ple2d.raw.plpar = dummy.plpar;
ar.ple2d.raw.par = dummy.par;

end

function q_trial = FindInsideJumps(eps)

% Check all raw 2d-profile points and save all points which satisfy that
% the chi2-value of the grid for the data point above and for the data
% point below is each lower than the raw profile point in between.
%
% These are potential outliers. Note that outliers in parameter
% direction are probably already dealt with with pleSmooth.

global ar

zq = ar.ple2d.smooth.zq;
q_trial = NaN(1,3);

for indy_tmp = 1:size(zq,1)
    
    %Get raw parameter profile ponts:
    xs = ar.ple2d.raw.plpar(:,indy_tmp);
    indxsraw = 1:length(xs);
    indxsraw = indxsraw(~isnan(xs));
    
    for indxraw_tmp = indxsraw
        
        chi2_tmp = ar.ple2d.raw.chi2(indxraw_tmp,indy_tmp); %Current chi2
        xs_tmp = sort([xs(indxraw_tmp),ar.ple2d.smooth.xq(1,:)]);
        
        % Find neighboring grid parameter indices:
        indx_down = find(xs_tmp == xs(indxraw_tmp),1)-1;
        indx_up = find(xs_tmp == xs(indxraw_tmp),1);
        if (indx_down == 0) || (indx_up == size(zq,2)+1)
            continue
        end
        % Find neighboring prediction indices:
        indy_down = indy_tmp - 1;
        indy_up = indy_tmp +1;
        if (indy_down == 0) || (indy_up == size(zq,1)+1)
            continue
        end
        
        % Get chi2 value for same parameter for data point below and data
        % point above. (Might be better to compare to raw profile values).
        if sum(isnan(zq(indy_down,[indx_down,indx_up]))) == 0
            chi2_down = interp1(ar.ple2d.smooth.xq(1,[indx_down,indx_up]),...
                zq(indy_down,[indx_down,indx_up]),xs(indxraw_tmp));
        else
            chi2_down = NaN;
        end
        if sum(isnan(zq(indy_up,[indx_down,indx_up]))) == 0
            chi2_up = interp1(ar.ple2d.smooth.xq(1,[indx_down,indx_up]),...
                zq(indy_up,[indx_down,indx_up]),xs(indxraw_tmp));
        else
            chi2_up = NaN;
        end
        
        % If the point in the middle exceeds both lower and upper point by
        % eps, it might be a candidate for smoothing.
        if (chi2_down + eps < chi2_tmp) && (chi2_up + eps < chi2_tmp)
            if isnan(q_trial(1,1))
                q_trial(1,1:3) = [xs(indxraw_tmp),indy_down,indy_up];
            else
                q_trial(end+1,1:3) = [xs(indxraw_tmp),indy_down,indy_up];
            end
        end
    end
end

% If no candidate has been found, return empty vector
if isnan(q_trial(1,1))
    q_trial = [];
end

end

