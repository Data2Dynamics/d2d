% Disturb model parameters
%
% arDisturb(strength)
%   strength:	pval+randn(0,strength)     [1]

function arDisturb(strength)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if(~exist('strength', 'var'))
    strength = 0.1;
end

ar.p(ar.qFit==1) = ar.p(ar.qFit==1) + randn(size(ar.p(ar.qFit==1))) * strength;
ar.p(ar.p<ar.lb) = ar.lb(ar.p<ar.lb);
ar.p(ar.p>ar.ub) = ar.ub(ar.p>ar.ub);

arChi2;