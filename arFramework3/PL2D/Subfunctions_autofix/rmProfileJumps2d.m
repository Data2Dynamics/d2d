function rmProfileJumps2d(findmode,indy,direc,lbub,indx)

% rmProfileJumps2d(findmode,lbub,indy,direc)
%
% This function fixes major discontinuities in the 2d-profile which occur
% between parameter profiles as a whole. This is done by recalculating
% the problematic profiles with new trial parameters (determined by a
% profile with neighboring prediction value). This is done iteratively until 
% no discontinuity is detected.
%
% indy:     Prediction index of profile from which starting parameter is
%               determined.
% direc:    Prediction direction in which profiles are fixed
% findmode: Determines how trial parameter is found. 
%               findmode = 1: Optimized for SmoothingPoints2d_Excess
%               findmode = 2: Optimized for SmoothingPoints2d_Bounds
%               findmode = 3: Optimized for SmoothingPoints2d_Inside
%               findmode = 4: Manual mode. Additionally specify smooth 
%                               parameter grid index where discontinuity occurs.
% lbub:     Use upper (1) or lower (-1) bound for reference. 
% indx:     Starting parameter grid index (only used for findmode = 4)
%
% Possible discontinuities are found by SmoothingPoints2d_Bounds,
% SmoothingPoints2d_Excess or SmoothingPoints2d_Inside. Note that they 
% share a set of same discontinuities which they find, but each method may 
% find unique discontinuities.
%
% See also: autofix2dprofile , SmoothingPoints2d_Excess
%               SmoothingPoints2d_Bounds , SmoothingPoints2d_Inside


global ar

try
    ar_old = arDeepCopy(ar);
    smooth2d;
    
    idpar = ar.ple2d.general.idpar;
    indy_tmp = indy;
    
    % The dummy struct contains the improved profiles and is saved into
    % ar.ple2d at the end of this function
    dummy.chi2 = ar.ple2d.raw.chi2;
    dummy.plpar = ar.ple2d.raw.plpar;
    dummy.par = ar.ple2d.raw.par;
    
    while ~((indy_tmp == 1) && direc  == -1) ...
            && ~((indy_tmp == size(ar.ple2d.raw.chi2,2)) && direc ==1)
        %% Find the trial parameter vector from the previous profile
        
        % This is done in order for smooth2d to be applicable.
        ar.ple2d.raw.chi2 = dummy.chi2;
        ar.ple2d.raw.plpar = dummy.plpar;
        ar.ple2d.raw.par = dummy.par;
        smooth2d;
        
        indy_prev = indy_tmp;
        indy_tmp = indy_tmp+direc;
        
        if exist('p_trial','var')
            p_trial_last = p_trial;
        else
            p_trial_last = [];
        end
        
        if findmode == 1
            % Use with SmoothingPoints2d_Excess
            %
            % Find starting profile parameter by starting from the minimum
            % of the previous profile excess:
            q_grid = ~isnan(ar.ple2d.smooth.zq);
            if lbub == -1
                indx_prev = find(q_grid(indy_prev,:),1,'first');
                indx_tmp = find(q_grid(indy_tmp,:),1,'first');
            elseif lbub == 1
                indx_prev = find(q_grid(indy_prev,:),1,'last');
                indx_tmp = find(q_grid(indy_tmp,:),1,'last');
            end
            if (indx_prev - indx_tmp == 0) || ...
                    ((indx_prev - indx_tmp < 0) && lbub == 1) || ...
                    ((indx_prev - indx_tmp > 0) && lbub == -1)
                break;
                % exit if the next profile exceeds the previous profile
            else
                indx_list = (indx_tmp+lbub):lbub:indx_prev;
                % For these parameter indices there is a profile excess
                chi2_tmp = ar.ple2d.smooth.zq(indy_prev,:);
                [~,subindx_tmp] = min(chi2_tmp(indx_list));
                indx_trial = indx_list(subindx_tmp);
            end
            
            exit_aux = ar.ple2d.config.autofix.excess_factor * ...
                icdf('chi2',ar.ple2d.config.bounds.alpha_bounds,1) ...
                            + min(ar.ple2d.smooth.zq(indy_tmp,:));  
            % exit_aux is needed to complete the mode-specifix exit
            % condition in fixprofile. Exit if starting profile value
            % exceeds the specified chi2-boundary.
        elseif findmode == 2
            % Use with SmoothingPoints2d_Bounds
            %
            % Parameter index is chosen such that the old 2d-profile has a
            % value at this point, which serves as a value to compare a new
            % trial optimum against.
            if lbub == -1
                indx_trial = find(~isnan(ar.ple2d.smooth.zq(indy_tmp,:)),1,'first');
            elseif lbub == 1
                indx_trial = find(~isnan(ar.ple2d.smooth.zq(indy_tmp,:)),1,'last');
            end            
            exit_aux = ar.ple2d.smooth.zq(indy_tmp,indx_trial) - ...
                            ar.ple2d.config.autofix.eps_outside;  
            % Exit if new optimum is not significantly better than
            % existing optimum
        elseif findmode == 3
            % Use with SmoothingPoints2d_Inside
            %
            % Parameter index is chosen such that the difference between
            % one profile and the next is maximal (likely discontinuity)
            [~,indx_trial] = max(ar.ple2d.smooth.zq(indy_tmp,:)...
                - ar.ple2d.smooth.zq(indy_prev,:));
            exit_aux = ar.ple2d.smooth.zq(indy_tmp,indx_trial) - ...
                            ar.ple2d.config.autofix.eps_inside;  
            % Exit if new optimum is not significantly better than
            % existing optimum                        
        elseif findmode == 4
            % Manual mode
            indx_trial = indx;
            exit_aux = ar.ple2d.smooth.zq(indy_tmp,indx_trial) - ...
                            ar.ple2d.config.autofix.eps_inside;                         
        else
            fprintf(['\n ERROR rmProfileJumps2d: Argument findmode is',...
                ' specified incorrectly. \n'])
            return
        end
        
        % Find whole corresponding parameter vector:
        x_trial = ar.ple2d.smooth.xq(1,indx_trial);
        p_trial = ptrialUpdate(dummy,indy_prev,indy_tmp,x_trial,p_trial_last);
        % The functions which find the discontinuities are built such that
        % p_trial can be updated this way at least for the first correction
        % step.
       
        %% Obtain new profile and merge with old one:
        
        % Calculate profile with new starting parameter:        
        q = fixprofile(indy_tmp,p_trial,exit_aux);
        
        % Save relevant results from calculation before resetting ar:
        ple1.chi2 = ar.ple2d.raw.chi2(:,indy_tmp);
        ple1.par = ar.ple2d.raw.par{indy_tmp};
        
        ple2.chi2 = (ar.ple.chi2s{idpar} ...
            - ar.ple2d.raw.chi2_min)';
        ple2.par = ar.ple.ps{idpar};
        
        % Reset ar after modification in fixprofile
        ar = ar_old;
        
        % exit conditions are defined by virtue of exit_aux and checked in 
        % fixprofile
        if q == 0
            break
        end
        
        % Merge all relevant fields of new profile with old one
        % Take the values of the profile which is lower than the other:
        [chi2,par] = mergeprofiles(idpar,ple1,ple2);
        
        dummy = StoreFixedResults(dummy,chi2,par,indy_tmp,idpar);
        
    end
catch exception
    fprintf(['ERROR rmProfileJumps2d: Resetting ar struct.',...
        ' Error message: \n %s \n Line: %s \n'] ,...
        exception.message, sprintf('%i, ',exception.stack.line));
    ar = ar_old;
end

ar.ple2d.raw.chi2 = dummy.chi2;
ar.ple2d.raw.plpar = dummy.plpar;
ar.ple2d.raw.par = dummy.par;

end



