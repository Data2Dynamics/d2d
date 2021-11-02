function q = fixprofile(indy,p_trial,exit_aux)
% This function calculates the profile with a new initial parameter, which
% is specified such that the profile starts in a better local optimum.

global ar

[m,d,idpred,iddata] = Set2dDataPlaceholder;
ar.model(m).data(d).yExp(iddata,idpred) = ar.ple2d.smooth.yq(indy,1);

% Initialize new profile starting point:
idpar = ar.ple2d.general.idpar;
ar.p = p_trial;
chi2 = UsePLEFit(idpar);

% Check whether new starting point has a lower chi^2 value than existing
% point
chi2_new = chi2 - ar.ple2d.raw.chi2_min;
del_new = chi2_new - min(ar.ple2d.smooth.zq(indy,:));
if (exit_aux - chi2_new > 0) 
    q = 1;
else
    % if custom exit condition is met, fix2dprofile terminates
    q = 0;
    return
end
 
% Calculate new profile:
SetPartialPLEInit(del_new,idpar);
ple(idpar);
pleSmooth;

% Can be used to compare old profile with new profile:
% hold on
% plot(ar.ple2d.raw.plpar(:,indy),ar.ple2d.raw.chi2(:,indy));
% plot(ar.ple.ps{idpar}(:,idpar),ar.ple.chi2s{idpar}...
%     -ar.ple2d.raw.chi2_min - ar.ple2d.smooth.offset);
% hold off

end



