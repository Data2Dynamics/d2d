% Select most parsimoneous model found by L1 regularization
% jks    relative parameters to be investigated by L1 regularization

function l1SelectOpt(jks)

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
ps = ar.L1ps;
ps_unpen = ar.L1ps_unpen;
chi2s_unpen = ar.L1chi2s_unpen;

parsgt0 = sum(abs(ps(:,jks)) > 1e-6,2);
parsgt0(1) = length(jks);

signifmat = nan(1,length(chi2s_unpen));
for j = 1:length(chi2s_unpen)
    signifmat(j) = chi2s_unpen(j) - chi2s_unpen(1) - icdf('chi2',.95,parsgt0(1)-parsgt0(j));
end
tmp = find(signifmat(1,:) < 0);
final_ind = tmp(end);

ar.p = ps_unpen(final_ind,:);
ar.type(jks) = 0;
ar.qFit(jks) = 1;
ar.qFit(jks(abs(ps(final_ind,jks)) <= 1e-6)) = 2;

ar.L1final_ind = final_ind;
ar.L1parsgt0 = parsgt0;
ar.L1signifmat = signifmat;

fprintf('Most parsimoneous model: %i / %i parameter(s) cell-type specific.\n',parsgt0(final_ind),length(jks))