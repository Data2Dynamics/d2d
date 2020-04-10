% p = arGetPars(pLabel, [qLog10])
%
% Get parameter value by matching label to ar.pLabel
% 
% pLabel	name of the parameter
% qLog10	logical to get log10 of parameter value [ar.qLog10]
%
% p         parameter value
%
% See also arGetParsPattern arPrint

function p = arGetPars(pLabel, qLog10)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~iscell(pLabel))
    pLabel = {pLabel};
end

qp = ismember(ar.pLabel, pLabel); %R2013a compatible

if(nargin < 2)
    qLog10 = ar.qLog10(qp);
end
if(length(qLog10)==1 && length(pLabel)>1)
    qLog10 = (zeros(size(pLabel)) + qLog10) == 1;
end

p = zeros(1,length(pLabel));
for j=1:length(pLabel)
    q = ismember(ar.pLabel, pLabel(j)); %R2013a compatible

    if(sum(q)==0)
        error('parameter %s not found', pLabel{j});
    end
    
    p(j) = ar.p(q);
    
    if(ar.qLog10(q) && ~qLog10(j))
        p(j) = 10^p(j);
    elseif(~ar.qLog10(q) && qLog10(j))
        p(j) = log10(p(j));
    end
end
