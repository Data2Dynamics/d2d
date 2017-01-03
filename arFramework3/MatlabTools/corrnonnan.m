function [RHO,PVAL] = corrnonnan(X, type)

if(~exist('type', 'var'))
    type = 'Pearson';
end

if(nargout<2 && sum(sum(isnan(X),2)==0)>0)
    RHO = corr(X, 'type', type);
    return
end

RHO = ones(size(X,2));
PVAL = ones(size(X,2));

for j1=1:size(X,2)
    for j2=(j1+1):size(X,2)
        x = X(:,[j1 j2]);
        x = x(sum(isnan(x),2)==0,:);
        if(size(x,1)>0)
            if(any(std(x)==0))
                if(j1==j2)
                    RHO(j1,j2) = 1;
                    PVAL(j1,j2) = 0;
                else
                    RHO(j1,j2) = 0;
                    RHO(j2,j1) = 0;
                    PVAL(j1,j2) = 1;
                    PVAL(j2,j1) = 1;
                end
            else
                [rho, pval] = corr(x, 'type', type);
                RHO(j1,j2) = rho(1,2);
                PVAL(j1,j2) = pval(1,2);
                RHO(j2,j1) = rho(2,1);
                PVAL(j2,j1) = pval(2,1);
            end
        else
            RHO(j1,j2) = nan;
            PVAL(j1,j2) = nan;
            RHO(j2,j1) = nan;
            PVAL(j2,j1) = nan;
        end
    end
end
