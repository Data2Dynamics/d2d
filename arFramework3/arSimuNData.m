% Simulate n data points for current parameter settings
%
% arSimuNData(npoints, nrep, logtimes, m, d)
%   npoints:    number of time points for simulation        [10]
%   nrep:       number of repetitions                       [3]
%   logtimes:   log-distribute time points                  [false]
%   m:          model index                    
%   d:          data index                    

function arSimuNData(npoints, nrep, logtimes, m, d)

global ar

if(isempty(ar))
    error('please initialize by arInit')
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

if(~exist('m','var'))
    for jm=1:length(ar.model)
        for jd=1:length(ar.model(jm).data)
            arSimuNData(npoints, nrep, logtimes, jm, jd)
        end
    end
    return
end
if(~exist('d','var'))
    for jd=1:length(ar.model(m).data)
        arSimuNData(npoints, nrep, logtimes, m, jd)
    end
    return
end

if(~logtimes)
    tpoints = transpose(linspace(ar.model(m).data(d).tLimExp(1), ar.model(m).data(d).tLimExp(2), npoints));
else
    if(ar.model(m).data(d).tLimExp(1) <= 0)
        tmin = 1e-3;
    else
        tmin = ar.model(m).data(d).tLimExp(1);
    end
    tpoints = transpose(logspace(log10(tmin), log10(ar.model(m).data(d).tLimExp(2)), npoints));
end

ntpoints = [];
for j=1:nrep
    ntpoints = [ntpoints; tpoints]; %#ok<AGROW>
end

arSimuData(ntpoints, m, d);
