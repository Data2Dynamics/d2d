% pleClearNaNs([jk]) 
% 
% Remove NaNs in ar.ple. Useful when calculation stopped and ple shall be extended
%
%   jk   [1:length(ar.ple.chi2s)]     remove NaN in ar.ple(jk)


function pleClearNaNs(jk)

global ar

if(isempty(ar.ple))
    error('PLE ERROR: please initialize')
end 

if(nargin<1)
    jk = 1:length(ar.ple.chi2s);
end

for i = jk
    ikeep = ~isnan(ar.ple.chi2s{i});
    if ~isempty(ikeep)
        ar.ple.samplesize(i) = sum(ikeep);
        ar.ple.chi2s{i} = ar.ple.chi2s{i}(ikeep);
        ar.ple.chi2sinit{i} = ar.ple.chi2sinit{i}(ikeep);
        ar.ple.chi2sviolations{i} = ar.ple.chi2sviolations{i}(ikeep);
        ar.ple.chi2spriors{i} = ar.ple.chi2spriors{i}(ikeep);
        ar.ple.chi2spriorsAll{i} = ar.ple.chi2spriorsAll{i}(ikeep);
        ar.ple.ps{i} = ar.ple.ps{i}(ikeep,:);
        ar.ple.psinit{i} = ar.ple.psinit{i}(ikeep,:);
        ar.ple.psinitstep{i} = ar.ple.psinitstep{i}(ikeep,:);
        ar.ple.gradient{i} = ar.ple.gradient{i}(ikeep,:);
        ar.ple.dpLast{i} = ar.ple.dpLast{i}(ikeep);
        ar.ple.dpStep{i} = ar.ple.dpStep{i}(ikeep,:);
    end
end