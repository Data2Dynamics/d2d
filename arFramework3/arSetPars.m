% arSetPars(pLabel, p, qFit, qLog10, lb, ub, type, meanp, stdp)
% 
% set parameter value and properties by label
%
% pLabel	name of the parameter
% 
%           Several parameter names can be provided as cell:
%           - either single values are provided and set for all parameters 
%           - or the length of the provided values has to coincide with the
%             number of parameter labels
% p			value of the  parameter(s)
% qFit		0=fixed, 1=fitted, 2=constant
% qLog10	0=normal, 1=log10 parameter values
% lb		lower parameter bound(s)
% ub		upper parameter bound(s)
% type		0=box prior, 1=normal prior
% meanp		mean of normal prior(s)
% stdp		standard deviation of normal prior(s)
% 
% 
% Alternative call:
% arSetPars(arStruct2)
% 
%   In this case, the respective argumens are taken from the provided
%   ar-struct.

function arSetPars(pLabel_or_ar, p, qFit, qLog10, lb, ub, type, meanp, stdp)

global ar

if(nargin==1 && isstruct(pLabel_or_ar))
    arIn    = pLabel_or_ar;
    pLabel  = arIn.pLabel;
    p       = arIn.p;
    qFit    = arIn.qFit;
    qLog10  = arIn.qLog10;
    lb      = arIn.lb;
    ub      = arIn.ub;
    type    = arIn.type;
    meanp   = arIn.mean;
    stdp    = arIn.std;    
else
    pLabel = pLabel_or_ar;
end
    

if(isempty(ar))
    error('please initialize by arInit')
end

if(~iscell(pLabel))
    pLabel = {pLabel};
end

if(nargin>1 && ~isempty(p) && length(p)==1)
    p = zeros(size(pLabel)) + p;
end
if(nargin>2 && ~isempty(qFit) && length(qFit)==1)
    qFit = zeros(size(pLabel)) + qFit;
end
if(nargin>3 && ~isempty(qLog10) && length(qLog10)==1)
    qLog10 = zeros(size(pLabel)) + qLog10;
end
if(nargin>4 && ~isempty(lb) && length(lb)==1)
    lb = zeros(size(pLabel)) + lb;
end
if(nargin>5 && ~isempty(ub) && length(ub)==1)
    ub = zeros(size(pLabel)) + ub;
end
if(nargin>6 && ~isempty(type) && length(type)==1)
    type = zeros(size(pLabel)) + type;
end
if(nargin>7 && ~isempty(meanp) && length(meanp)==1)
    meanp = zeros(size(pLabel)) + meanp;
end
if(nargin>8 && ~isempty(stdp) && length(stdp)==1)
    stdp = zeros(size(pLabel)) + stdp;
end

for j=1:length(pLabel)
    q = ismember(ar.pLabel, pLabel(j)); %R2013a compatible
    if(sum(q)==1)
        if(nargin>1 && ~isempty(p))
            ar.p(q) = p(j);
            if(nargin>4 && ~isempty(lb))
                if(p(j) < lb(j))
                    warning('trying to set p < lb   [%.4f < %.4f (%d: %s)]', p(j), lb(j), j, pLabel{j});
                    ar.p(q) = lb(j);
                end
            else
                if(p(j) < ar.lb(q))
                    warning('trying to set p < lb   [%.4f < %.4f (%d: %s)]', p(j), ar.lb(q), find(q), pLabel{j});
                    ar.p(q) = ar.lb(q);
                end
            end
            if(nargin>5 && ~isempty(ub))
                if(p(j) > ub(j))
                    warning('trying to set p > ub   [%.4f > %.4f (%d: %s)]', p(j), ub(j), j, pLabel{j});
                    ar.p(q) = ub(j);
                end
            else
                if(p(j) > ar.ub(q))
                    warning('trying to set p > ub   [%.4f > %.4f (%d: %s)]', p(j), ar.ub(q), find(q), pLabel{j});
                    ar.p(q) = ar.ub(q);
                end
            end
        end
        if(nargin>2 && ~isempty(qFit))
            ar.qFit(q) = qFit(j);
        end
        if(nargin>3 && ~isempty(qLog10))
            ar.qLog10(q) = qLog10(j);
        end
        if(nargin>4 && ~isempty(lb))
            ar.lb(q) = lb(j);
        end
        if(nargin>5 && ~isempty(ub))
            ar.ub(q) = ub(j);
        end
        if(nargin>6 && ~isempty(type))
            ar.type(q) = type(j);
        end
        if(nargin>7 && ~isempty(meanp))
            ar.mean(q) = meanp(j);
        end
        if(nargin>8 && ~isempty(stdp))
            ar.std(q) = stdp(j);
        end
    elseif(sum(q)==0)
        arFprintf(1, 'arSetPars: parameter %s not found\n', pLabel{j});
    else
        error('multiple parameters %s!?', pLabel{j});
    end
end
