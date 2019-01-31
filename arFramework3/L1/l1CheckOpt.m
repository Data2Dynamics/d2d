%% DEPRECATED
% l1CheckOpt(jks)
% 
% L1 check optimum
% After most parsimoneous model is found, do single PLE for each non-specific
% parameter included in the fitted set to check if others cross 0 at
% re-optimization. If so, this parameter could also be cell-type specific.
% 
% jks             [ar.L1jks]
%                 indices of the fold-factor parameters to be investigated by L1
%                 regularization 

function l1CheckOpt(jks)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
        error('please initialize by l1Init, run l1Scan, and l1Unpen')
    end
end
jks = ar.L1jks;

fixed_jks = find(abs(ar.L1ps(ar.L1final_ind,jks)) <= 1e-4);

ar.p = ar.L1ps_unpen(ar.L1final_ind,:);
ar.type(jks) = 0;
ar.qFit(jks) = 1;
ar.qFit(jks(fixed_jks)) = 2;

arPLEInit
for jk = 1:length(jks)
    ar.ple.p = ar.L1ps_unpen(ar.L1final_ind,:);
    ar.ple.q_fit(jks) = 1;
    ar.ple.q_fit(jks(fixed_jks)) = 0;
    ar.ple.q_fit(jks(jk)) = 1;
    ar.qFit(jks) = 1;
    ar.qFit(jks(fixed_jks)) = 2;
    ar.qFit(jks(jk)) = 1;
    
%     arFit(true)
%     arPLEInit
    
    do_plotting = ar.ple.showCalculation;
    ar.ple.showCalculation = false;
    ple(jks(jk),50,0.1,0.1,0.1)
    ar.ple.showCalculation = do_plotting;
    
    ar.L1psPLE{jks(jk)} = ar.ple.ps{jks(jk)};
    ar.L1chi2sPLE{jks(jk)} = ar.ple.chi2s{jks(jk)};
    not_profiled = setdiff(jks,jks(jk));
    pleSigns = sign(ar.L1psPLE{jks(jk)}(:,not_profiled));
    exchange_jk = find(max(pleSigns)-min(pleSigns) == 2);
    if ~isempty(exchange_jk)
        fprintf('Parameter #%i: ''%s'' could be exchanged by #%i: ''%s''\n',jks(jk),ar.pLabel{jks(jk)},not_profiled(exchange_jk),ar.pLabel{not_profiled(exchange_jk)});
    end
end