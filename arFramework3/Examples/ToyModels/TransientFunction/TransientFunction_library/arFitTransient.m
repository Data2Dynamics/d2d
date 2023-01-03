% arFitTransient
% 
% The major function for fitting the transient function. It calls arFit and
% tries upward- and downward responses.
% 
% This function fits serveral times:
%   - upward direction
%   - downward direction
% ar.fit_transient.doReduction determines whether backward elimination is
% applied:
%   - for LRT for testing toffset=0
%   - for LRT for testing f = constant
%   - for LRT for testing offset=0
% 
%   The signums are first set to +1. Then one of the other signum is set
%   to -1 and it is checked whether the fit improves.
% 
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   !  Since the bounds depend on the signum, this function also   !
%   !  temporarily changes upper and lower bounds.                             !
%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function arFitTransient
global ar

if isfield(ar.fit_transient,'doReduction')
    doReduction = ar.fit_transient.doReduction;
else
    doReduction = true;
end

ar.lb = ar.fit_transient.bounds.lb;
ar.ub = ar.fit_transient.bounds.ub;

indsig = ar.fit_transient.indp.signum;

if sum(ar.p(ar.qFit==1)>ar.ub(ar.qFit==1))>0
    indfit = find(ar.qFit==1);    
    arPrint(indfit(find(ar.p(indfit)>ar.ub(indfit))))
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
fits.up = arCopyFit(ar); % all signum positive 
arCalcMerit
chi2best = arGetMerit;

% Now try fits which each signum negative.
% Since the functions are independent, one does not need all combis.
for i=1:length(indsig) 
    ind_arp = ar.fit_transient.(ar.pLabel{indsig(i)}).ind_arp;
    if sum(ar.qFit(ind_arp)==1) > 0 % some paraemters of this transient function are fitted
        lbstart = ar.lb+0.0;
        ubstart = ar.ub+0.0;
        
        pstart = ar.p+0.0;
        
        ar.p = p0;
        ar.lb(indsig(i)) = ar.fit_transient.boundsNeg.lb(indsig(i));
        ar.ub(indsig(i)) = ar.fit_transient.boundsNeg.ub(indsig(i));
        ar.p(indsig(i)) = -1;
        arFit
        arCalcMerit;
        fits.down{i} = arCopyFit( ar);  
        
        if(chi2best < arGetMerit) % no improvement
            fprintf('=> %i th signum: regulation in upward direction was better. Reset p.\n',i)
            ar.p = pstart; % reset ar.p
            ar.lb = lbstart; % reset lb
            ar.ub = ubstart; % reset ub
            arCalcMerit
        else
            fprintf('=> %i th signum: regulation in downward direction. Keep new p.\n',i)            
            % take better fit and update ar:
            chi2best = arGetMerit;
        end
    else
        fprintf('arFitTransient: Parameter %s is fixed: SKIPPED.\n',ar.pLabel{indsig(i)});
    end
end
 
fits.doReduction = doReduction;
if doReduction
    %% Try setting toffset to zero (LRT with current best fit, although the number of data points is arbitrary and the error is fitted)
    last.p = ar.p+0.0;
    last.qFit = ar.qFit+0.0;
    last.qLog10 = ar.qLog10+0.0;
    
    ar.p(ar.fit_transient.indp.toffset) = ar.lb(ar.fit_transient.indp.toffset);
    ar.qFit(ar.fit_transient.indp.toffset) = 0;
    ar.qLog10(ar.fit_transient.indp.toffset) = 0;
    
    arFit
    fits.toffsetZero = arCopyFit(ar);
    if(chi2best < arGetMerit-icdf('chi2',.99,length(ar.fit_transient.indp.toffset))) % no improvement
        ar.qLog10 = last.qLog10;
        ar.qFit = last.qFit;
        ar.p = last.p;
        fprintf('LRT rejected: toffset ~= lb\n');
        arCalcMerit
    else
        chi2best = arGetMerit;
        fprintf('LRT: toffset = lb\n');
    end
    %% Try using only a constant
    last.p = ar.p+0.0;
    last.qFit = ar.qFit+0.0;
    last.qLog10 = ar.qLog10+0.0;
    
    ar.qFit(ar.fit_transient.indp.toffset) = 0;
    ar.qFit(ar.fit_transient.indp.amp_sust) = 0;
    ar.qFit(ar.fit_transient.indp.amp_trans) = 0;
    ar.qFit(ar.fit_transient.indp.timescale_sust) = 0;
    ar.qFit(ar.fit_transient.indp.timescale_trans) = 0;
    
    ar.p(ar.fit_transient.indp.toffset) = ar.lb(ar.fit_transient.indp.toffset);
    ar.p(ar.fit_transient.indp.amp_sust) = ar.lb(ar.fit_transient.indp.amp_sust);
    ar.p(ar.fit_transient.indp.amp_trans) = ar.lb(ar.fit_transient.indp.amp_trans);
    ar.p(ar.fit_transient.indp.timescale_sust) = ar.lb(ar.fit_transient.indp.timescale_sust);
    ar.p(ar.fit_transient.indp.timescale_trans) = ar.lb(ar.fit_transient.indp.timescale_trans);
    
    
    arFit
    fits.offsetZero = arCopyFit(ar);
    pTest = [ar.fit_transient.indp.toffset,...
        ar.fit_transient.indp.amp_sust,...
        ar.fit_transient.indp.amp_trans,...
        ar.fit_transient.indp.timescale_sust,...
        ar.fit_transient.indp.timescale_trans];
    pTest = intersect(pTest,find(last.qFit==1));
    np = length(pTest);
    
    if(chi2best < arGetMerit-icdf('chi2',.99,np)) % no improvement
        ar.p = last.p;
        ar.qFit = last.qFit;
        ar.qLog10 = last.qLog10;
        fprintf('LRT rejected: Trajectore NOT constant.\n');
    else
        chi2best = arGetMerit;
        fprintf('LRT: Trajectory constant => Use only offset parameter.\n');
    end
    
    
    %% Try setting offset to zero
    last.p = ar.p+0.0;
    last.qFit = ar.qFit+0.0;
    last.qLog10 = ar.qLog10+0.0;
    
    ar.p(ar.fit_transient.indp.offset) = 0;
    ar.qFit(ar.fit_transient.indp.offset) = 0;
    ar.qLog10(ar.fit_transient.indp.offset) = 0;
    
    arFit
    fits.offsetZero = arCopyFit(ar);
    if(chi2best < arGetMerit-icdf('chi2',.99,length(ar.fit_transient.indp.offset))) % no improvement
        ar.qLog10 = last.qLog10;
        ar.qFit = last.qFit;
        ar.p = last.p;
        fprintf('LRT rejected: offset ~= 0\n');
    else
        chi2best = arGetMerit;
        fprintf('LRT: offset = 0\n');
    end
end
%%

ar.fit_transient.fits = fits;
