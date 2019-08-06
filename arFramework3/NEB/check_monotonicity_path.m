function [is_path_OK] = ct_check_mono_path(chi2s, tolerance_chi2)

[~,xx] = min(chi2s);

if xx > 1 && xx < length(chi2s)
    is_path_OK = issorted(chi2s(1:xx),'descend') && issorted(chi2s(xx:end),'ascend');
else
    is_path_OK = issorted(chi2s,'descend')|| issorted(chi2s,'ascend');
end


if is_path_OK == 0 
    % further  Analysis
    
    %left hand side
    if issorted(chi2s(1:xx),'descend') == 0 
    %issue is left
        
        % first deviation from monotony left
        wi = 0;
        n = xx-1;
        while issorted(chi2s(n:xx),'descend')==1
            wi = wi +1;
            n = xx - wi;
        end
        
        % from xx_left to xx(minimum) all ok, thus check only
        % 1:xx_left
        xx_left = n+1;   
         
        if max(chi2s(1:xx_left)) -  min(chi2s(1:xx_left)) < tolerance_chi2
            % is flat left, within tolerance
            is_path_OK_left = 1;
        else % not flat within tolerance
            %
            [is_path_OK_left, corrected_path_l]  = ...
                check_stepwise_tolerance_left(chi2s(1:xx_left), tolerance_chi2);
        end
        
     else % links ist OK
        is_path_OK_left = 1;

    end %left
    
    %right hand side
    if issorted(chi2s(xx:end),'ascend') == 0 
    %issue is right
        
        % first deviatiom from monotony right
        wi = 0;
        n = xx+1;
        while issorted(chi2s(xx:n),'ascend')==1
            wi = wi + 1;
            n = xx + wi;
        end
        
        % from xx(minimum) to xx_right all ok, thus check only
        % xx_right:end
        xx_right = n;   
        
        if max(chi2s(xx_right:end)) -  min(chi2s(xx_right:end)) < tolerance_chi2
            % is flat right, within tolerance
            is_path_OK_right = 1;
        else % not flat within tolerance
            % try interpolate 
            
            [is_path_OK_right, corrected_path_r]  = ...
                check_stepwise_tolerance_right(chi2s(xx_right-1:end), tolerance_chi2);

        end
        
    else % right hand side is OK
        is_path_OK_right = 1;
        
    end % right
    
    %check if everything is OK now
    if is_path_OK_right == 1 && is_path_OK_left == 1
        is_path_OK = 1;
    end
end 

end
