% result = arIdentifiablityTest_WhichP
%
% arIdentifiablityTest_WhichP checks, which parameters have to be reoptimized if
% parameter p_name is fixed to p_val. An L1-penalty is set to ar.p
% 
% The current implementation only checks
% ar.IdentifiabilityTest.p_nonIdentifiable{1}. This should/will be
% generalized in more detail.
% 
%   Example:
% arLoad                    % e.g. Raia
% arFit
% arIdentifiablityTest
% arIdentifiablityTest_WhichP
% 
% See also arIdentifiablityTest

function result = arIdentifiablityTest_WhichP
global ar

result = struct;

if ~isfield(ar,'IdentifiabilityTest')
    warning('arIdentifiablityTest_WhichP.m: ... aborted. Perform arIdentifiabilityTest first!')
else
    
    [~,p_index] = intersect(ar.pLabel,ar.IdentifiabilityTest.p_nonIdentifiable{1});
    p_index = p_index(1);
    fprintf('%s is non-identifiable. Looking for related parameters ...\n',ar.pLabel{p_index});
    
    itype0 = find(ar.type==0 & ar.qFit==1);
    not_testable = setdiff(find(ar.qFit==1),itype0);
    if ~isempty(not_testable)
        fprintf('No testable due to penalties:\n');
        fprintf('  ''%s''\n',ar.pLabel{not_testable});
    end
    fprintf('\n');
    i_L1 = setdiff(itype0,p_index);
    
    mean_old = ar.mean;
    std_old = ar.std;
    
    if isempty(i_L1)
        fprintf('No parameters without penalties.\n');
    else
        
        ar.type(i_L1) = 3;
        ar.mean(i_L1) = ar.p(i_L1);
        sd = 1*ar.IdentifiabilityTest.radius;
        ar.std(i_L1) = sd;
        
        qfit0 = ar.qFit;
        p0 = ar.p;
        
        ar.qFit(p_index) = 0;
        ar.p(p_index) = ar.IdentifiabilityTest.p(p_index);
        arFit(true);
        dp = ar.p - p0;
        
        result.std = ar.std;
        result.mean = ar.mean;
        result.type = ar.type;
        result.i_L1 = i_L1;
        result.p_index = p_index;
        result.p = ar.p;
        result.dp = dp;
        
        ar.p = p0;
        ar.qFit = qfit0;
        ar.mean = mean_old;
        ar.std  = std_old;
        ar.type(i_L1) = 0;
        
        index_related = i_L1(abs(dp(i_L1)) > (0.01+0.01*sd));
        result.index_related = index_related;
        result.pLabel_related = ar.pLabel(index_related);
        
        if ~isempty(index_related)
            fprintf('Parameters in the non-identifiability-manifold of ''%s'':\n',ar.pLabel{p_index})
            fprintf('  ''%s''\n',ar.pLabel{index_related});
        else
            fprintf('Parameter ''%s'' is the only parameter in the non-identifiability-manifold.\n',ar.pLabel{p_index})
        end
    end
    
end

