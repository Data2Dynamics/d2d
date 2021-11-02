function SetPartialPLEInit(del_new,idpar)

% SetPartialPLEInit(del_new,idpar)
%
% This function sets the profile threshold such that profile is only 
% calculated as far as it needs to be sampled, since the default is a fixed 
% chi2-distance from the starting point (which is not necessarily the minimum).

global ar

arPLEInit(true,true,2,0)
ar.ple.showCalculation = 0;

ar.ple.dchi2_point = icdf('chi2',ar.ple2d.config.bounds.alpha_bounds,1);
% Since ple couples the stepsize to the threshold, the relative step size
% is increased if threshold becomes too small.
if del_new > 0.7*ar.ple.dchi2_point  
        ar.ple.relchi2stepincrease(idpar) = 2*ar.ple2d.config.gen2d.relchi;
end

% Since the actual threshold which is used for ple calculations is also
% controlled by ar.ple.dchi2_point, it needs to modified with an absolute
% inccrease if the new threshold would become too small.
if del_new > 0.8*ar.ple.dchi2_point
    ar.ple.dchi2_point = abs(ar.ple.dchi2_point-del_new) + 0.2*ar.ple.dchi2_point;
elseif del_new > 0.6*ar.ple.dchi2_point 
    ar.ple.dchi2_point = abs(ar.ple.dchi2_point-del_new) + 0.15*ar.ple.dchi2_point;
elseif del_new > 0.4*ar.ple.dchi2_point 
    ar.ple.dchi2_point = abs(ar.ple.dchi2_point-del_new) + 0.1*ar.ple.dchi2_point;
elseif del_new > 0.2*ar.ple.dchi2_point 
    ar.ple.dchi2_point = abs(ar.ple.dchi2_point-del_new) + 0.05*ar.ple.dchi2_point;
end

end

