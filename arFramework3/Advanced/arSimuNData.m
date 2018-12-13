% arSimuNData(npoints], [nrep], [logtimes], [m], [d], [randomseed])
% 
% Simulate n data points for each observable for current parameter
% settings by calling arSimuData.
%
%   npoints:    number of time points for simulation        [10]
%   nrep:       number of repetitions                       [3]
%   logtimes:   log-distribute time points                  [false]
%   m:          model index                    
%   d:          data index              
%   randomseed  random seed for noise generation
% 
% The following commands are used to select time-points:
% linspace(ar.model(m).data(d).tLimExp(1), ar.model(m).data(d).tLimExp(2), npoints) % logtimes = true
% logspace(log10(tmin), log10(ar.model(m).data(d).tLimExp(2)), npoints) % logtimes = true
% 
% This function is useful when creating realistic designs.

function arSimuNData(npoints, nrep, logtimes, m, d, randomseed)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if(~exist('npoints','var'))
    npoints = 10;
end

if(~exist('randomseed','var'))
    randomseed = [];
end

if(~exist('nrep','var'))
    nrep = 3;
end
if(~exist('logtimes','var'))
    logtimes = false;
end

if(~exist('m','var') || isempty(m))
    for jm=1:length(ar.model)
        for jd=1:length(ar.model(jm).data)
            arSimuNData(npoints, nrep, logtimes, jm, jd, randomseed)
        end
    end
    return
end
if(~exist('d','var') || isempty(d))
    for jd=1:length(ar.model(m).data)
        arSimuNData(npoints, nrep, logtimes, m, jd, randomseed)
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

arSimuData(m, d, ntpoints, randomseed);
