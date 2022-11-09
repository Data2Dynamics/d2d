function [Sigma, dSigma_dp] = arCalcCovarianceMatrix(idm,idd,idy)

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
    sdCov = 0;
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

varianceMatrix = sd*sd';

deltatMatrix = abs(repmat(tExp,1,length(tExp))-repmat(tExp',length(tExp),1));

Sigma = varianceMatrix .* sdCov.^deltatMatrix;

if nargout > 1
    dSigma_dp = NaN([size(Sigma),length(ar.model(idm).data(idd).p)]);
    systdExpSimu = ar.model(idm).data(idd).systdExpSimu(~idnan,:,:);
    systdExpSimu_trafo = arTrafoParameters(systdExpSimu,idm,idd,1);
    for idp = 1:length(ar.model(idm).data(idd).p)
        ssd = systdExpSimu_trafo(:,idy,idp);

        idpCov = find(strcmp(ar.model(idm).data(idd).p,ar.model(idm).data(idd).pcov(idy) ));

        if idp==idpCov
            if ar.qLog10(ar.model(idm).data(idd).pcovLink(idy,:)==1)
                dSigma_dp(:,:,idp) = varianceMatrix.*deltatMatrix.*sdCov.^(deltatMatrix-1)*sdCov*log(10); % dSigma_dsdCov
            else
                dSigma_dp(:,:,idp) = varianceMatrix.*deltatMatrix.*sdCov.^(deltatMatrix-1); % dSigma_dsdCov
            end
        else
            dSigma_dp(:,:,idp) = (ssd * sd' + sd * ssd') .* sdCov.^deltatMatrix;
        end
    end
end