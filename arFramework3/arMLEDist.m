% check distribution of MLE
%
% arMLEDist(ntimes, npoints, nrep, logtimes)
%
%   ntimes      number of realizations                      [100]
%   npoints:    number of time points for simulation        [10]
%   nrep:       number of repetitions                       [3]
%   logtimes:   log-distribute time points                  [false]

function arMLEDist(ntimes, npoints, nrep, logtimes)

global ar
global mledist

if(~exist('ntimes','var'))
    ntimes = 1000;
end
if(~exist('npoints','var'))
    npoints = 10;
end
if(~exist('nrep','var'))
    nrep = 3;
end
if(~exist('logtimes','var'))
    logtimes = false;
end

if(~isfield(ar, 'pTrue'))
    ar.pTrue = ar.p;
end

mledist.ntimes = ntimes;
mledist.npoints = npoints;
mledist.nrep = nrep;

mledist.chi2 = nan(1,ntimes);
mledist.chi2err = nan(1,ntimes);
mledist.chi2fit = nan(1,ntimes);
mledist.pdiff = nan(ntimes, length(ar.p));

h = waitbar(0, 'Please wait...');

for j=1:ntimes
    if(mod(j,10)==0)
        h = waitbar(j/ntimes, h, 'Please wait...');
    end
    
    arReset(true);
    arSimuData(npoints, nrep, logtimes);
    try
        arFit(true);
        
        mledist.chi2(j) = ar.chi2 + 0;
        mledist.chi2err(j) = ar.chi2err + 0;
        mledist.chi2fit(j) = ar.chi2fit + 0;
        mledist.pdiff(j,:) = ar.p - ar.pTrue + 0;
    catch expception
        fprintf('#%i ERROR: %s\n', j, expception.message);
    end
end

close(h);

mledist.ndata = ar.ndata;
mledist.npara = sum(sum(abs(ar.sres),1)>0 & ar.qFit==1);
mledist.fiterrors = ar.config.fiterrors;

arReset(true);
