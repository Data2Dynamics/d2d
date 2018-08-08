% This function fits serveral times.
% 
%   The signums are first set to +1. Then one ofter the other signum is set
%   to -1 and it is checked whether the fit improves.
% 
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   !  Since the founds depend on the signum, this function also   !
%   !  changes upper and lower bounds.                             !
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function arFitTransient
global ar

ar.lb = ar.fit_transient.bounds.lb;
ar.ub = ar.fit_transient.bounds.ub;

indsig = ar.fit_transient.indp.signum;

if sum(ar.p(ar.qFit==1)>ar.ub(ar.qFit==1))>0
    error('ar.p>ar.ub')
end
if sum(ar.p(ar.qFit==1)<ar.lb(ar.qFit==1))>0
    error('ar.p<ar.lb')
end
if sum(ar.ub(ar.qFit==1)<ar.lb(ar.qFit==1))>0
    error('ar.ub<ar.lb')
end

p0 = ar.p + 0.0;

ar.p(indsig) = 1; % set all signums to one
arFit
fits.up = ar; % all signum positive
arCalcMerit
chi2best = arGetMerit;

% Now try fits which each signum negative.
% Since the functions are independent, one does not need all combis.
for i=1:length(indsig) 
    ind_arp = ar.fit_transient.(ar.pLabel{indsig(i)}).ind_arp;
    if sum(ar.qFit(ind_arp)==1) > 0 % some paraemters of this transient function are fitted
        lbstart = ar.lb;
        ubstart = ar.ub;
        
        pstart = ar.p+0.0;
        
        ar.p = p0;
        ar.lb(indsig(i)) = ar.fit_transient.boundsNeg.lb(indsig(i));
        ar.ub(indsig(i)) = ar.fit_transient.boundsNeg.ub(indsig(i));
        ar.p(indsig(i)) = -1;
        arFit
        fits.down{i} = ar;
        
        if(chi2best < arGetMerit) % no improvement
            fprintf('=> %i th signum: regulation in upward direction was better. Reset p.\n',i)
            ar.p = pstart; % reset ar.p
            ar.lb = lbstart; % reset lb
            ar.ub = ubstart; % reset ub
        else
            fprintf('=> %i th signum: regulation in downward direction. Keep new p.\n',i)
            
            % take better fit and update ar:
            arCalcMerit;
            chi2best = arGetMerit;
        end
    else
        fprintf('arFitTransient: Parameter %s is fixed: SKIPPED.\n',ar.pLabel{indsig(i)});
    end
end

%% Try setting toffset to zero (LRT with current best fit, although the number of data points is arbitrary and the error is fitted)
last.p = ar.p;
last.qFit = ar.qFit;
last.qLog10 = ar.qLog10;

ar.p(ar.fit_transient.indp.toffset) = ar.lb(ar.fit_transient.indp.toffset);
ar.qFit(ar.fit_transient.indp.toffset) = 0;
ar.qLog10(ar.fit_transient.indp.toffset) = 0;

arFit
fits.toffsetZero = ar;
ar.qFit = last.qFit;
ar.qLog10 = last.qLog10;
if(chi2best < arGetMerit-icdf('chi2',.99,length(ar.fit_transient.indp.toffset))) % no improvement
    ar.p = last.p;
    fprintf('LRT rejected: toffset ~= lb\n');
else
    fprintf('LRT: toffset = lb\n');
end
%% Try setting offset to zero 
last.p = ar.p;
last.qFit = ar.qFit;
last.qLog10 = ar.qLog10;

ar.p(ar.fit_transient.indp.offset) = 0;
ar.qFit(ar.fit_transient.indp.offset) = 0;
ar.qLog10(ar.fit_transient.indp.offset) = 0;

arFit
fits.offsetZero = ar;
ar.qFit = last.qFit;
ar.qLog10 = last.qLog10;
if(chi2best < arGetMerit-icdf('chi2',.99,length(ar.fit_transient.indp.offset))) % no improvement
    ar.p = last.p;
    fprintf('LRT rejected: offset ~= 0\n');
else
    fprintf('LRT: offset = 0\n');
end
%%

ar.fit_transient.fits = fits;
