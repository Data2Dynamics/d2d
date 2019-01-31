% grplasSelectOpt([jks], [method])
% 
% Select most parsimoneous model found by group-lasso, i.e. the
% chosen penalization strength.
% 
% jks       [ar.L1jks]
%           indices of the fold-factor parameters to be investigated by L1
%           regularization 
% method    ['LRT'] is default
%           'BIC' is a possible alternative
%           The method used for definining the most parsimonious model
%           
% See also l1SelectOpt

function grplasSelectOpt(jks,method)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~isfield(ar,'grplas'))
    error('please initialize by grplasInit')
end

if(~exist('jks','var') || isempty(jks))
    if(~isfield(ar.grplas,'jks') || isempty(ar.grplas.jks))
        error('please initialize by grplasInit, run grplasScan, and grplasUnpen')
    end
end

if(~exist('method','var')|| isempty(method))
    method = 'LRT';
end


jks = ar.grplas.jks;
ps = ar.grplas.ps;
ps_unpen = ar.grplas.ps_unpen;
chi2s_unpen = ar.grplas.chi2s_unpen;
chi2s_lam0 = ar.grplas.lam0chi2s;

parsgt0 = sum(abs(ps(:,jks)) > ar.grplas.thresh,2);
parsgt0 = parsgt0';
parsgt0lam0 = length(jks);

signifmat = [];

if strcmpi(method,'LRT')
    
    for j = 1:length(chi2s_unpen)
        signifmat(j) = chi2s_unpen(j) - chi2s_lam0 - icdf('chi2',.95,parsgt0lam0-parsgt0(j));
    end

    tmp = find(signifmat(1,:) < 0 | parsgt0 == parsgt0lam0);
    if ~isempty(tmp)
        final_ind = tmp(end);
    else
        final_ind = 1;
    end
elseif strcmpi(method,'BIC')
    estpars = sum(ar.type ~= 5) + parsgt0;
    signifmat = log(ar.ndata) * estpars + chi2s_unpen;
    [~, final_ind]  = min(signifmat);
end
    
ar.p = ps_unpen(final_ind,:);
ar.type(jks) = 0;
ar.qFit(jks) = 1;
ar.qFit(jks(abs(ps(final_ind,jks)) <= ar.grplas.thresh)) = 2;
%ar.p(jks(abs(ps(final_ind,jks)) <= ar.grplas.thresh)) = 0;

ar.grplas.final_ind = final_ind;
ar.grplas.parsgt0 = parsgt0;
ar.grplas.signifmat = signifmat;
ar.grplas.seltype = method;

fprintf('Most parsimoneous model: %i / %i parameter(s) cell-type specific.\n',parsgt0(final_ind),length(jks))