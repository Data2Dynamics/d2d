function a = meannonnan(b, dim)
if(isvector(b))
    q = ~isnan(b);
    a = mean(b(q));
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(b,1),1);
        for j = 1:size(b,1)
            a(j) = meannonnan(b(j,:));
        end
    else
        a = nan(1,size(b,2));
        for j = 1:size(b,2)
            a(j) = meannonnan(b(:,j));
        end
    end
end