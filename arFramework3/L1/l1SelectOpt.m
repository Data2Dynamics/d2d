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

% Calculate differences between fc parameters
if ar.L1DiffPen_activate
    diffIds = ar.L1DiffPen_diffs;
    ps2 = bsxfun(@minus, ps(:,diffIds(:,1)), ps(:,diffIds(:,2)));
    jks2 = 1:size(diffIds,1);
    excl2 = cell(length(chi2s_unpen),1);
end

parsgt0lam0 = length(jks);

signifmat = nan(1,length(chi2s_unpen));

if strcmpi(method,'LRT')
    for j = 1:length(chi2s_unpen)
        if ar.L1DiffPen_activate == true
            %excl = jks(abs(ps(j,jks)) <= ar.L1thresh); % find fcs that are zero
            %excl2{j} = jks2(abs(ps2(j,jks2)) <= ar.L1thresh); % find fc diffs that are zero
            %exclFromExcl2 = diffIds(excl2{j},:); % find fcs that belong to zero fc diffs
            %reduceParsGt0By = sum(sum(ismember(exclFromExcl2,excl),2) == 0);
            
            reduceParsGt0By = 0;
            redDiff = false(size(diffIds,1),1);
            for k = 1:size(diffIds,1)
                id1 = diffIds(k,1);
                id2 = diffIds(k,2);
                if abs(ps2(j,k)) < ar.L1thresh
                    if (abs(ps(j,id1)) > ar.L1thresh && abs(ps(j,id2)) > ar.L1thresh)
                        
                        rows1 = any(diffIds(redDiff,:) == id1, 2);         % Find pairs for which parsgt0 is reduced that contain id1
                        rows2 = any(diffIds(redDiff,:) == id2, 2);
                        partners1 = diffIds(rows1, diffIds(rows1) ~= id1); % Find the partner in the pair that is not id1
                        partners2 = diffIds(rows2, diffIds(rows2) ~= id2);

                        if sum(ismember(partners1,partners2)) == 0         % If any partners match, skip reduction of parsgt0 
                            redDiff(k) = 1;
                            reduceParsGt0By = reduceParsGt0By + 1;
                        end
                    end
                end
            end
            
            % redDiffIds = diffIds(reducedDiff,:);
            % id = redDiffIds == id1
            % test(any(id,2),id(any(id,2),:))
            
            if reduceParsGt0By ~= 0
                fprintf('l1SelectOpt: lambda #%i: Reduced parsgt0 by %i\n',j,reduceParsGt0By)
            end
            parsgt0(j) = parsgt0(j) - reduceParsGt0By; %%% CUSTOM
        end
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

if ar.L1DiffPen_activate 
    fprintf('The following delta is zero: %i\n', excl2{final_ind})
end

ar.p = ps_unpen(final_ind,:);
ar.type(jks) = 0;
ar.qFit(jks) = 1;
ar.qFit(jks(abs(ps(final_ind,jks)) <= ar.L1thresh)) = 2;
if ar.L1DiffPen_activate
    diffZeros = ar.L1DiffPen_constrSwitches(jks2(abs(ps2(final_ind,:)) <= ar.L1thresh));
    ar.qFit(diffZeros) = 2;
    ar.p(diffZeros) = 1;
    ar.qLog10(diffZeros) = 0;
end
% ar.p(jks(abs(ps(final_ind,jks)) <= ar.L1thresh)) = 0;

ar.L1final_ind = final_ind;
ar.L1parsgt0 = parsgt0;
ar.L1signifmat = signifmat;
ar.L1seltype = method;

fprintf('Most parsimoneous model: %i / %i parameter(s) cell-type specific.\n',parsgt0(final_ind),length(jks))