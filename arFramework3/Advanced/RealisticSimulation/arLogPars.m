function arLogPars

global ar

id = ar.p>0 & ar.lb>0 & ~ar.qLog10;
idb0 = ar.p>0 & ar.lb==0 & ~ar.qLog10;
idp0 = find(ar.p==0 & ~ar.qLog10);

ar.p(id|idb0) = log10(ar.p(id|idb0));
ar.lb(id) = log10(ar.lb(id));
ar.ub(id|idb0) = log10(ar.ub(id|idb0));
ar.qLog10(id|idb0) = 1;

try
    arSimu(false,true);
end
if length(idp0)>0
    low = nan(length(idp0),1);
    for i=1:length(idp0)
        idx = contains(ar.model.x,ar.pLabel{idp0(i)}(6:end));
        if sum(idx)==1
	mi = min(log10(ar.model.condition.xFineSimu(ar.model.condition.xFineSimu(:,idx)>0,idx)));
	if isempty(mi)
	   mi = -5;
	end
            low(i) = mi;
        elseif sum(idx)>1
            mi = min(min(log10(ar.model.condition.xFineSimu(ar.model.condition.xFineSimu>0)))); 
	if isempty(mi)
	   mi = -5;
	end
            low(i) = mi;
        else
	low(i) = -5;
        end
    end

    ar.p(idp0) = low;
    ar.lb(idp0) = low -2;
    ar.ub(idp0) = low + 2;
    ar.qLog10(idp0) = 1;
    
    ar.lb(idb0) = min(low);
else
    ar.lb(idb0) = min(min(log10(ar.model.condition.xFineSimu(ar.model.condition.xFineSimu>0))));
end


if any(ar.p(idb0)<ar.lb(idb0))
    ar.lb(idb0) = ar.p(idb0)-2;
    warning(['arLogPars: Lower bound of ' ar.pLabel{idb0} ' set to ' num2str(ar.p(idb0)-2) '.'])
end
if ~all(id|idb0)
    warning('%s ',['arLogPars: Following parameters were not logged due to negative parameter or negative boundary value: ' ar.pLabel{~(id|idb0)} ' '])
end