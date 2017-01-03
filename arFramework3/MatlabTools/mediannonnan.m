function a = mediannonnan(b, dim)
if(isvector(b))
    q = ~isnan(b);
    a = median(b(q));
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(b,1),1);
        for j = 1:size(b,1)
            a(j) = mediannonnan(b(j,:));
        end
    else
        a = nan(1,size(b,2));
        for j = 1:size(b,2)
            a(j) = mediannonnan(b(:,j));
        end
    end
end