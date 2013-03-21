% get parameter value by label
%
% p = arGetPars(pLabel, qLog10)
% 
% pLabel	name of the parameter
% qLog10	0=normal, 1=log10 parameter values

function p = arGetPars(pLabel, qLog10)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~iscell(pLabel))
    pLabel = {pLabel};
end

qp = ismember(ar.pLabel, pLabel);

if(nargin < 2)
    qLog10 = ar.qLog10(qp);
end
if(length(qLog10)==1 && length(pLabel)>1)
    qLog10 = (zeros(size(pLabel)) + qLog10) == 1;
end

p = zeros(1,length(pLabel));
for j=1:length(pLabel)
    q = ismember(ar.pLabel, pLabel(j));

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
