function a = wstd(c, dim)
if(isvector(c) && ~exist('dim','var'))
    a = 1/sum(c.^-2);
else
    if(exist('dim','var') && dim == 2)
        a = nan(size(c,1),1);
        for j = 1:size(c,1)
            a(j) = meannonnan(c(j,:));
        end
    else
        a = nan(1,size(c,2));
        for j = 1:size(c,2)
            a(j) = meannonnan(c(:,j));
        end
    end
end