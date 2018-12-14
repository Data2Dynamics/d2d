% arFitEvA2([optimizer_index, silent])
%
% user EvA2 global optimizers
%
%     optimizer_index: [4] 
%     1: STD_ES         Standard (15,50) Evolutionary Strategy
%     2: CMA-ES         with Covariance Matrix Adaption
%     3: STD_GA         Standard Genetic Algorithm
%     4: PSO            Particle Swarm Optimization with constriction
%     5: DE             Differential Evolution
%     6: TRIBES         adaptive PSO
%     7: RANDOM         Random Monte Carlo Search
%     8: HILLCL         Hill-Climbing 
%     9: CBN_ES         Cluster-based niching ES 
%     10: CHILL         Clustering Hill-Climbing 
%     11: IPOP-CMA-ES 
%     12: CBN_GA        Cluster-based niching GA 
%     13: PBIL
%
%     silent - boolean to suppress output [false]

function arFitEvA2(optimizer_index, silent)

global ar

if(~exist('optimizer_index','var'))
    optimizer_index = 4; % starts a default PSO run 
end
if(nargin<2)
    silent = false;
end

ar.fevals = 0;
ar.ps_errors = [];

JI = JEInterface(@testfun, 'double', [ar.lb(ar.qFit==1);ar.ub(ar.qFit==1)], ...
    [ar.lb(ar.qFit==1);ar.ub(ar.qFit==1)]);%, [], ...
%     'Display', 'off');
% 	'Display', ar.config.optim.Display);
% 	  'TolFun', ar.config.optim.TolFun, ...
%     'TolX', ar.config.optim.TolX, ...

% disp(getDesc(JI, optimizer_index));

JI = optimize(JI, optimizer_index); 
sol = getResult(JI);
ar.p(ar.qFit==1) = sol;

ar.fit.message = getMessage(JI);

if(~silent)
    arCalcMerit;
else
    arCalcMerit(false);
end



function z = testfun(x, ~)
global ar

try
    arCalcMerit(false, x);
catch id
    disp(id.message);
    ar.ps_errors(end+1,:) = ar.p;
end
z = arGetMerit('chi2fit');