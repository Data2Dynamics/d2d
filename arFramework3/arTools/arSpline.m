function [y, dydp] = arSpline(xp,yp,x)

y = spline(xp,yp,x);

if(nargout>1)
    dydp = nan(length(x),length(yp));
    for j=1:length(yp)
        yptmp = zeros(size(yp));
        yptmp(j) = 1;
        dydp(:,j) = spline(xp,yptmp,x);
    end
end