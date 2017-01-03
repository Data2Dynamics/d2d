% Calculate sensitivities for fitting by finite difference approximations 
% of residuals ar.res in ar.sresFD and of constraints ar.constr in ar.sconstrFD
%
% arFiniteDifferences(dp)    
%   dp:     parameter variation             [1e-6]


function arFiniteDifferencesGlobal(dp)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('dp', 'var'))
    dp = 1e-6;
end

pRef = ar.p;
arChi2(false,ar.p(ar.qFit==1));

if(~isempty(ar.res))
    ar.sresFD = ar.res' * ones(1,length(pRef)) + 0;
else
    ar.sresFD = [];
end

if(~isempty(ar.constr))
    ar.sconstrFD = ar.constr' * ones(1,length(pRef)) + 0;
else
    ar.sconstrFD = [];
end

arWaitbar(0);
for jp=1:length(ar.pLabel)
    arWaitbar(jp, length(ar.pLabel));
    
    % disturb p(jp)
    ar.p = pRef;
    ar.p(jp) = ar.p(jp) + dp;
    arChi2(false,ar.p(ar.qFit==1));
    
    if(~isempty(ar.res))
        ar.sresFD(:,jp) = (ar.res' - ar.sresFD(:,jp)) / dp + 0;
    end
    
    if(~isempty(ar.constr))
        ar.sconstrFD(:,jp) = (ar.constr' - ar.sconstrFD(:,jp)) / dp + 0;
    end
end

ar.p = pRef;
arChi2(true,ar.p(ar.qFit==1));

arWaitbar(-1);

