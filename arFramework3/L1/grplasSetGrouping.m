%GRPLASSETGROUPING Assigns a grouping to ar
%   Detailed explanation goes here
% 
% grplasSetGrouping([grouping], [jks], weights )
% 
% grouping      Definition of the parameters groups which are
%               zero/non-zero simultaneously:
%               [ones(size(jks))] is default and means all jks are in one group
%               Alternatively a cell array of parameter names can be used.
%               Alternative call: 'alt...' (not yet documented).
%               Alternative call: 'group...' (not yet documented).
% 
% jks           indices of the fold-factor parameters to be investigated by L1
%               regularization 
%               [find(ar.type == 5)] is default 
% 
% weights       Scalar, vector or matrix defining group weights
%               [1] 

function grplasSetGrouping( grouping, jks, weights )
global ar

if(~exist('jks','var') || isempty(jks))
    indices = 1:length(ar.p);
    jks = indices(ar.type == 5);
end

if(~exist('grouping','var') || isempty(grouping))
    grouping = ones(size(jks));
elseif iscell(grouping)
    tmpgrp = nan(size(jks));
    for i = 1:length(jks)
        for j = 1:length(grouping)
            if contains(ar.pLabel(jks(i)),grouping{j})
                tmpgrp(i) = j;
                break
            end
        end
        if isnan(tmpgrp(i))
            pname = ar.pLabel(jks(i));
            s = sprintf('Cannot assign group to %s',pname{1});
            error(s)
        end
    end
    grouping = tmpgrp;
elseif strcmp( grouping(1:3), 'alt' )
    [n, b] = str2num( grouping(4:end) );
    if b
        
        jind = 1:length(jks);
        grouping = mod(jind-1, n)+1;
        
    else
        error('Input not recognized, please use alt2, alt12 etc.')
        
    end
elseif strcmp( grouping(1:5) , 'group' )
    [n, b] = str2num( grouping(6:end) );
    if b
        i = 1;
        j = 1;
        while (i < length(jks))
            tmpgrp(i:(i +n-1)) = j;
            i = i + n;
            j = j + 1;
        end
        grouping = tmpgrp;
    else
        error('Input not recognized, please use group2, group3 etc.')
    end
elseif strcmp(grouping , 'single' )
    warning('please use l1 methods')
    grouping = 1:length(jks);
end

        

if(~exist('weights','var') || isempty(weights))
    weights = 1;
end



ar.grplas.grouping = zeros(size(ar.p));
ar.grplas.grouping(jks) = grouping;
ar.grplas.groups = unique(grouping);
ar.grplas.Ngroups = length(ar.grplas.groups);
ar.grplas.groupinds = {};
for g = ar.grplas.groups
    ar.grplas.groupinds{g} = find(ar.grplas.grouping == g);
end

ar.grplas.A = double(diag(ar.type==5));

if isscalar(weights)
    % the same weights for all groups
    ar.grplas.A = ar.grplas.A * weights;
elseif (isvector(weights))
    % one weight per group
    if isequal(size(weights), size(jks))
        % one weight per parameter
        ar.grplas.A(jks,jks) = diag(abs(weights));
    elseif isequal(size(weights),size(ar.grplas.groups))
        tmpA = zeros(1,length(ar.type));
        for i = 1:length(ar.grplas.groups)
            gind = (ar.type == 5) & (ar.grplas.grouping == ar.grplas.groups(i));
            tmpA(gind) = abs(weights(i));
        end
        ar.grplas.A = diag(tmpA);
    end
elseif ismatrix(weights)
    
    try
        chol((weights+weights')./2);
    catch
        error('Please provide a positive-definit matrix as weights')
    end
    
    if isequal(size(weights),[length(jks) length(jks)])
        ar.grplas.A(jks,jks) = weights;
    else
        for g = 1:ar.grplas.Ngroups
            Ng = length(ar.grplas.groupinds{g});
            if isequal(size(weights),[Ng Ng])
                ar.grplas.A(ar.grplas.groupinds{g},ar.grplas.groupinds{g}) ...
                    = weights;
            else
                ar.grplas.A(ar.grplas.groupinds{g},ar.grplas.groupinds{g}) ...
                    = eye(Ng);
            end
        end
    end
    
end

ar.grplas.A = sparse(ar.grplas.A);
% ar.grplas.weights = ones(1,ar.grplas.Ngroups) * weights;

ar.grplas.eig = eig((ar.grplas.A+ar.grplas.A')./2);
    

end

