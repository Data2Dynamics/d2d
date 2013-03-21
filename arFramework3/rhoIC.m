% Statistics for model selection

function rhoIC
global ar


nFreePar = sum(ar.qFit==1);

AIC = ar.chi2fit + ar.ndata*log(2*pi) + 2*nFreePar;
BIC = ar.chi2fit + ar.ndata*log(2*pi) + nFreePar*log(ar.ndata);
AICc = ar.chi2fit + ar.ndata*log(2*pi) + 2*nFreePar*(nFreePar+1)/(ar.ndata-nFreePar-1);

fprintf('AIC = %f, AICc = %f, BIC = %f\n',AIC, AICc, BIC);