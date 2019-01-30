% arSample([N], [parindex], [mode], [range])
% 
% Samples the likelihood in parameter space in directions of parindex.
% 
%   N           number of sampling steps within the range       [100]
%   parindex    vector of indices of parameters to sample. 
%               Maximal length is three.                        [1]
%   mode        1 for ub to lb                                  [1]
%               2 for range of PLE
%               3 for user defined range
%   range       give the range where to sample (mode = 3): 
%               range is a cell array of vectors with length two indicating 
%               lower and upper bound per parameter in parindex. [{[] [] []}]
% 
% Results are stored in ar.sampling.
% For three parameters (parindex) and N = 100, this leads to 10^6 merit function
% evaluations.
% 
% See also arPlotChi2s arScan arScanChi2s

function arSample(N, parindex, mode, range)

global ar

if(~exist('N','var'))
    N = 100;
end
if(~exist('parindex','var'))
    parindex = 1;
end
if(~exist('mode','var'))
    mode = 1;
end
if(mode == 3 && ~exist('range','var'))
    error('range argument required');
elseif(mode ~= 3)
    range = {[] [] []};
end

ar.sampling.index = parindex;

pTrue = ar.p;
arWaitbar(0);
if(length(parindex)==1)
    ar.sampling.ps = {makerange(N,parindex, mode, range{1})};
    ar.sampling.chi2s = nan(1,N);
    for j=1:N
        arWaitbar(j, N);
        ar.p(parindex) = ar.sampling.ps{1}(j);
        try %#ok<TRYNC>
            arCalcMerit(false);
        end 
        ar.sampling.chi2s(j) = arGetMerit('chi2fit');
    end
elseif(length(parindex)==2)
    ar.sampling.ps = {makerange(N,parindex(1), mode, range{1}) makerange(N,parindex(2), mode, range{2})};
    ar.sampling.chi2s = nan(N);
    for j1=1:N
        for j2=1:N
            arWaitbar(j2+(j1-1)*N, N*N);
            ar.p(parindex) = [ar.sampling.ps{1}(j1) ar.sampling.ps{2}(j2)];
            try %#ok<TRYNC>
                arCalcMerit(false);
            end
            ar.sampling.chi2s(j2,j1) = arGetMerit('chi2fit');
        end
    end
    dx = ar.sampling.ps{1}(2)-ar.sampling.ps{1}(1);
    dy = ar.sampling.ps{2}(2)-ar.sampling.ps{2}(1);
    ar.sampling.marginalized{1} = sum(exp(-0.5*ar.sampling.chi2s),1)*dy;
    ar.sampling.marginalized{2} = sum(exp(-0.5*ar.sampling.chi2s),2)*dx;
elseif(length(parindex)==3)
    ar.sampling.chi2s = nan(N,N,N);
else
    close(h);
    error('maximum of three dimensions allowed');
end
ar.p = pTrue;
arCalcMerit(false);
arWaitbar(-1);


function p = makerange(N,parindex, mode, range)
global ar

if(mode==1)
    p = linspace(ar.lb(parindex), ar.ub(parindex), N);
elseif(mode==2)
    p = linspace(min(ar.ple.ps{parindex}(:,parindex)), ...
        max(ar.ple.ps{parindex}(:,parindex)), N);
elseif(mode==3)
    p = linspace(range(1), range(2), N);
end
