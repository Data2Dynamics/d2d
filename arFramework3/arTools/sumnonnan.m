function a = sumnonnan(b, dim)
if(isvector(b))
    q = ~isnan(b);
    a = sum(b(q));
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(b,1),1);
        for j = 1:size(b,1)
            a(j) = sumnonnan(b(j,:));
        end
    else
        a = nan(1,size(b,2));
        for j = 1:size(b,2)
            a(j) = sumnonnan(b(:,j));
        end
    end
end