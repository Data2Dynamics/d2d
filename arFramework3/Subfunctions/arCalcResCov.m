%function [] = arCalcResCov(sensi)
% This function calculates the transformed residuals used for the
% autocorrelation error model from the original residuals.
%
%       sensi   should sensitivities sres, sreserr be calculated?
%               [sensi=1]
%
% The following fields in the ar struct are filled by this function:
%      ar.model(...).data(...).resCov
%      ar.model(...).data(...).reserrCov
%      ar.model(...).data(...).chi2cov
%      ar.model(...).data(...).chi2err_cov
%      ar.model(...).data(...).sresCov
%      ar.model(...).data(...).sreserrCov
%
% see also arCollectResCov, arFit, arGetMerit

function [] = arCalcResCov(sensi)
if ~exist('sensi','var') || isempty(sensi)
    sensi = 1;
end

global ar

for idm = 1:length(ar.model)
    if(isfield(ar.model(idm), 'data'))
        for idd = 1:length(ar.model.data)
            
            checkTimePoints(idm,idd);
            
            [ar.model(idm).data(idd).resCov, ar.model(idm).data(idd).chi2cov] = ...
                calcResCov(idm,idd);
            [ar.model(idm).data(idd).reserrCov, ar.model(idm).data(idd).chi2err_cov] = ...
                calcReserrCov(idm,idd);
            if sensi
                ar.model(idm).data(idd).sresCov = ...
                    calcSensiResCov(idm,idd);
                ar.model(idm).data(idd).sreserrCov = ...
                    calcSensiReserrCov(idm,idd);
            end
        end
    end
end

end


function [resCov,chi2cov] = calcResCov(idm,idd)

global ar

res = ar.model(idm).data(idd).res;% nTimes x nObs
resCov = nan(size(res));
for idy = 1:size(res,2)
    res_idy = res(:,idy);
    res_idy_minusOne = [0; res_idy(1:end-1)];
    phis = getPhis(idm,idd,idy);
    
    resCov(:,idy) = 1./sqrt(1-phis.^2) .* (res_idy - phis .* res_idy_minusOne);
end
chi2cov = sum(resCov.^2,1);

end

function [reserrCov, chi2err_cov] = calcReserrCov(idm,idd)

global ar

reserr = ar.model(idm).data(idd).reserr;% nTimes x nObs
reserrCov = nan(size(reserr));
chi2err_cov = ar.model(idm).data(idd).chi2err;
for idy = 1:size(reserr,2)
    reserr_idy = reserr(:,idy);
    phis = getPhis(idm,idd,idy);
    
    reserrCov_idy_2 = reserr_idy.^2 + log(1-phis.^2);
    if(sum(reserrCov_idy_2(:) < 0)>0)
        error('arCalcResCov/calcReserrCov: covariance error residual get imaginary. Increase ar.config.add_c!');
    else 
        reserrCov(:,idy) = sqrt(reserrCov_idy_2); %reserr_idy includes ar.config.add_c already
        chi2err_cov(idy) = chi2err_cov(idy) + sum(log(1-phis.^2),1);
    end
end

end

function [sresCov] = calcSensiResCov(idm,idd)

global ar

sres = ar.model(idm).data(idd).sres;% nTimes x nObs x nPars
sresIdxMinusOne = nan(size(sres));
sresIdxMinusOne(2:end,:,:) = sres(1:end-1,:,:);
sresIdxMinusOne(1,:,:) = 0;

res = ar.model(idm).data(idd).res;% nTimes x nObs

sresCov = nan(size(sres));
for idy = 1:size(res,2)
    res_idy = res(:,idy);
    res_idy_minusOne = [0; res_idy(1:end-1)];
    [phis, dphis_dcov] = getPhis(idm,idd,idy);
    
    idpCov = find(strcmp(ar.model(idm).data(idd).p,ar.model(idm).data(idd).pcov(idy)));
    for idp = 1:length(ar.model(idm).data(idd).p)
        if idp==idpCov
            sresCov(:,idy,idp) = ...
                dphis_dcov .* ( ...
                res_idy .* phis .* (1-phis.^2).^(-1.5) - ...
                res_idy_minusOne .* (1-phis.^2).^(-1.5) ...
                );
        else
            sresCov(:,idy,idp) = ...
                1./sqrt(1-phis.^2) .* (sres(:,idy,idp) - phis .* sresIdxMinusOne(:,idy,idp));
        end
    end
end

end

function [sreserrCov] = calcSensiReserrCov(idm,idd)

global ar

reserr = ar.model(idm).data(idd).reserr;% nTimes x nObs
reserrCov = ar.model(idm).data(idd).reserrCov;% nTimes x nObs
sreserr = ar.model(idm).data(idd).sreserr;% nTimes x nObs x nPars

sreserrCov = (reserr./reserrCov) .* sreserr;
for idy = 1:size(sreserr,2)
    idpCov = find(strcmp(ar.model(idm).data(idd).p,ar.model(idm).data(idd).pcov(idy)));
    if ~isempty(idpCov)
        [phis, dphis_dcov] = getPhis(idm,idd,idy);
        sreserrCov(:,idy,idpCov) = sreserrCov(:,idy,idpCov) - ...
            1./reserrCov(:,idy) .* phis./(1-phis.^2) .* dphis_dcov;
    end
end

end

function [phis,dphis_dcov] = getPhis(idm,idd,idy)

global ar

sd = ar.model(idm).data(idd).ystdExpSimu(:,idy);
tExp = ar.model(idm).data(idd).tExp;

idnan = isnan(sd);
sd(idnan) = [];
tExp(idnan) = [];

if any(ar.model(idm).data(idd).pcovLink(idy,:))
    if ar.qLog10(ar.model(idm).data(idd).pcovLink(idy,:)==1)
        sdCov = 10.^ar.p(ar.model(idm).data(idd).pcovLink(idy,:)==1);
    else
        sdCov = ar.p(ar.model(idm).data(idd).pcovLink(idy,:)==1);
    end
else
    phis = zeros(length(tExp),1);
    dphis_dcov = zeros(length(tExp),1);
    return
end

if sdCov >= 1
    warning('Covariance parameter %s has to be smaller than 1 on the linear scale.', ar.pLabel{ar.model(idm).data(idd).pcovLink(idy,:)==1})
    if ar.qLog10(ar.model(idm).data(idd).pcovLink(idy,:)==1)
        arSetPars(ar.pLabel{ar.model(idm).data(idd).pcovLink(idy,:)==1}, log10(0.99), [], [], [], log10(0.99))
    else
        arSetPars(ar.pLabel{ar.model(idm).data(idd).pcovLink(idy,:)==1}, 0.99, [], [], [], 0.99)
    end
    sdCov = 0.99;
    disp('Parameter value and upper bound was set to 0.99 on linear scale.')
end

phis = [0; sdCov.^abs(diff(tExp))];
if nargout==2
    if ar.qLog10(ar.model(idm).data(idd).pcovLink(idy,:)==1)
        dphis_dcov = log(10) .* [0; abs(diff(tExp)) .* sdCov.^(abs(diff(tExp)))];
    else
        dphis_dcov = [0; abs(diff(tExp)) .* sdCov.^(abs(diff(tExp))-1)];
    end
end

end


function checkTimePoints(idm,idd)

global ar

if length(unique(ar.model(1).data(1).tExp)) ~= length(ar.model(1).data(1).tExp)
    error('arCalcResCov: covariance error model not applicable for duplicate time points!')
end

end
