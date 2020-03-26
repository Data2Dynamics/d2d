function pleIterativeUntilConvergence(deltaBetter,nRestartMax,varargin)
if(~exist('deltaBetter','var') || isempty(deltaBetter))
    deltaBetter = 0.02;
end
if(~exist('nRestartMax','var') || isempty(nRestartMax))
    nRestartMax = 20;
end


nRestart = 0;
arCalcMerit;
chi2old = arGetMerit;

arPLEInit
ple(varargin{:})
ckBestPFromPle2;
arCalcMerit;
chi2 = arGetMerit;
arSave('current')

while chi2<chi2old-deltaBetter && nRestart<nRestartMax
    chi2old = chi2;
    nRestart = nRestart+1;
    fprintf('\npleIntativeUntilConvergence.m: %ith restart (out of %i)\n\n',nRestart,nRestartMax);

    arFit
    arPLEInit
    ple(varargin{:})
    ckBestPFromPle2;
    arCalcMerit;
    chi2 = arGetMerit;
    arSave('current')
end


