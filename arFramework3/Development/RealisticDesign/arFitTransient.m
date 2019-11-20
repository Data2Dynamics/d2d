function arFitTransient(n)

global ar
% 
[bounds,boundsNeg] = DefaultLbUbTransient;

% in.lb = ar.lb;
% in.ub = ar.ub;

ar.model.data.logfitting = 0;
ar.model.data.logplotting = 0;

ar.lb = bounds.lb;
ar.ub = bounds.ub;
ar.p(ar.fit_transient.indp.signum) = 1;
arFitLHS(n)
fit.up = ar;

ar.lb = boundsNeg.lb;
ar.ub = boundsNeg.ub;
ar.p(ar.fit_transient.indp.signum) = -1;
arFitLHS(n)
fit.down = ar;

if(fit.up.chi2fit<fit.down.chi2fit)
    better = 'up';
    disp('=> Regulation in upward direction.')
else
    better = 'down';    
    disp('=> Regulation in downward direction.')
end

% take better fit and update ar:
update_fields = {'chi2fit','chi2','chi2err','chi2prior','p'};
for i=1:length(update_fields)
    ar.(update_fields{i}) = fit.(better).(update_fields{i});
end
arSimu(false, true);

% ar.lb = in.lb;
% ar.ub = in.ub;
 
