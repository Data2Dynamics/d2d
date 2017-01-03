function a = stdnonnan(b, dim)
if(isvector(b))
    q = ~isnan(b);
    a = std(b(q));
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(b,1),1);
        for j = 1:size(b,1)
            a(j) = stdnonnan(b(j,:));
        end
    else
        a = nan(1,size(b,2));
        for j = 1:size(b,2)
            a(j) = stdnonnan(b(:,j));
        end
    end
end