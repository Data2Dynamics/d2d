function a = zscorenonnan(b, dim)
if(isvector(b))
    q = ~isnan(b);
    a = nan(size(b));
    a(q) = zscore(b(q));
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(b));
        for j = 1:size(b,1)
            a(j,:) = zscorenonnan(b(j,:));
        end
    else
        a = nan(size(b));
        for j = 1:size(b,2)
            a(:,j) = zscorenonnan(b(:,j));
        end
    end
end