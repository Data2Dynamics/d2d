function y = arSplineFit(xd, yd, xp, x, initialSlope)

if(~exist('x','var'))
    x = [];
end
if(~exist('initialSlope','var'))
    initialSlope = [];
end

yp0 = zeros(size(xp));
lb = zeros(size(xp))-5;
ub = zeros(size(xp))+5;
opts = optimoptions('lsqnonlin','Display','off');
yp = lsqnonlin(@(yp) mymerrit(xp, yp, xd, yd, initialSlope), yp0, lb, ub, opts);

if(isempty(x))
    y = spline(xp,yp);
else
    y = spline(xp,yp,x);
end



function F = mymerrit(xp, yp, xd, yd, initialSlope)

pp = spline(xp,yp);
yp = ppval(pp,xd); 
F = yp(:)-yd(:);

if(~isempty(initialSlope))
    ppdt = arSplineDer(pp);
    F(end+1) = ppval(ppdt,xp(1)) - initialSlope;
end