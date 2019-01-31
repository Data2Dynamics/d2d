% l1SelectOpt([jks], [method])
% 
% Select most parsimoneous model found by L1 regularization, i.e. the
% chosen penalization strength.
% 
% jks       [ar.L1jks]
%           indices of the fold-factor parameters to be investigated by L1
%           regularization 
% method    ['LRT'] is default
%           'BIC' is a possible alternative
%           The method used for definining the most parsimonious model
%           

function l1SelectOpt(jks,method)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar,'L1jks') || isempty(ar.L1jks))
        error('please initialize by l1Init, run l1Scan, and l1Unpen')
    end
end

if(~exist('method','var')|| isempty(method))
    method = 'LRT';
end

jks = ar.L1jks;
ps = ar.L1ps;
ps_unpen = ar.L1ps_unpen;
chi2s_unpen = ar.L1chi2s_unpen;
chi2s_lam0 = ar.L1lam0chi2s;

parsgt0 = sum(abs(ps(:,jks))> ar.L1thresh,2);
parsgt0lam0 = length(jks);

signifmat = nan(1,length(chi2s_unpen));

if strcmpi(method,'LRT')    
    for j = 1:length(chi2s_unpen)
        signifmat(j) = chi2s_unpen(j) - chi2s_lam0 - icdf('chi2',.95,parsgt0lam0-parsgt0(j));
    end

    tmp = find(signifmat(1,:) < 0 | parsgt0' == parsgt0lam0);
    if ~isempty(tmp)
        final_ind = tmp(end);
    else
        final_ind = 1;
    end
elseif strcmpi(method,'BIC')
    estpars = sum(ar.type ~= 5) + parsgt0';
    signifmat = log(ar.ndata) * estpars + chi2s_unpen;
    [~, final_ind]  = min(signifmat);
else
    error('method %s is not implemented.',method);
end
    

ar.p = ps_unpen(final_ind,:);
ar.type(jks) = 0;
ar.qFit(jks) = 1;
ar.qFit(jks(abs(ps(final_ind,jks)) <= ar.L1thresh)) = 2;
% ar.p(jks(abs(ps(final_ind,jks)) <= ar.L1thresh)) = 0;

ar.L1final_ind = final_ind;
ar.L1parsgt0 = parsgt0;
ar.L1signifmat = signifmat;
ar.L1seltype = method;

fprintf('Most parsimoneous model: %i / %i parameter(s) cell-type specific.\n',parsgt0(final_ind),length(jks))