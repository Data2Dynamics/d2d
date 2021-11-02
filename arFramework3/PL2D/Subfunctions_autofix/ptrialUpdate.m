function p_trial = ptrialUpdate(dummy,indy_prev,indy,x_trial,p_trial_last)

% Find parameter vector of the previous profile at the right place

try
    xs_tmp = dummy.plpar(:,indy_prev);
    ps_tmp = dummy.par{indy_prev};
    ps_tmp = ps_tmp(~isnan(xs_tmp),:);
    xs_tmp = xs_tmp(~isnan(xs_tmp));
    p_trial = interp1(xs_tmp,ps_tmp,x_trial);
catch
    fprintf('\n p_trial could not be updated for indy = %i \n',indy)
    p_trial = p_trial_last;
end
if sum(isnan(p_trial)) > 0
    fprintf('\n p_trial could not be updated for indy = %i \n',indy)
    p_trial = p_trial_last;
end

end

