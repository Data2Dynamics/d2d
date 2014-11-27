% Remove NaNs in pleGlobals
% Useful when calculation stopped and ple shall be extended
%
% pleClearNaNs([i])
%
% i:                    i'th parameter, see pwInfo

function pleClearNaNs(jk)

global pleGlobals;

if(isempty(pleGlobals))
    error('PLE ERROR: please initialize')
end 

if(nargin<1)
    jk = 1:length(pleGlobals.chi2s);
end

for i = jk
    ikeep = ~isnan(pleGlobals.chi2s{i});
    if ~isempty(ikeep)
        pleGlobals.samplesize(i) = sum(ikeep);
        pleGlobals.chi2s{i} = pleGlobals.chi2s{i}(ikeep);
        pleGlobals.chi2sinit{i} = pleGlobals.chi2sinit{i}(ikeep);
        pleGlobals.chi2sviolations{i} = pleGlobals.chi2sviolations{i}(ikeep);
        pleGlobals.chi2spriors{i} = pleGlobals.chi2spriors{i}(ikeep);
        pleGlobals.ps{i} = pleGlobals.ps{i}(ikeep,:);
        pleGlobals.psinit{i} = pleGlobals.psinit{i}(ikeep,:);
        pleGlobals.psinitstep{i} = pleGlobals.psinitstep{i}(ikeep,:);
        pleGlobals.gradient{i} = pleGlobals.gradient{i}(ikeep,:);
    end
end