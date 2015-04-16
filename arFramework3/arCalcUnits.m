function arCalcUnits(m)

global ar

if(~exist('m','var'))
    m = 1;
end

for jx = 1:length(ar.model(m).x)
    
end

return

ar.model(m).pUnits = cell(size(ar.model(m).p));
for jp = 1:length(ar.model(m).p)
    %% volumes
    q = ismember(ar.model(m).pc, ar.model(m).p{jp});
    if(sum(q) == 1)
        ar.model(m).pUnits(jp) = ar.model(m).c(q);
        continue
    elseif(sum(q) > 1)
        error('double matches');
    end
    
    %% initial conditions
    q = ismember(ar.model(m).px0, ar.model(m).p{jp});
    if(sum(q) == 1)
        ar.model(m).pUnits(jp) = ar.model(m).x(q);
        continue
    elseif(sum(q) > 1)
        error('double matches');
    end
    
    %% remaining dynamic parameters
    q = ismember(ar.model(m).pv, ar.model(m).p{jp});
    if(sum(q) == 1)
        q = ismembern(ar.model(m).pvs,ar.model(m).p{jp});
        vs = ar.model(m).fv(q);
        
        ptmp = '';
        
        ar.model(m).pUnits{jp} = ptmp;
        continue
    elseif(sum(q) > 1)
        error('double matches');
    end
end


function q = ismembern(c,l)
q = logical(size(c));
for j=1:length(c)
    q(j) = sum(ismember(c{j},l))>0;
end
    

