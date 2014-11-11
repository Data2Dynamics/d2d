function y = arSplineFit(xd, yd, xp, x, initialSlope, curveConstr)

if(~exist('x','var'))
    x = [];
end
if(~exist('initialSlope','var'))
    initialSlope = [];
end
if(~exist('curveConstr','var'))
    curveConstr = [];
end

yp0 = zeros(size(xp));
lb = zeros(size(xp))-5;
ub = zeros(size(xp))+5;
opts = optimoptions('lsqnonlin','Display','off');
yp = lsqnonlin(@(yp) mymerrit(xp, yp, xd, yd, initialSlope, curveConstr), yp0, lb, ub, opts);

if(isempty(x))
    y = spline(xp,yp);
else
    y = spline(xp,yp,x);
end



function F = mymerrit(xp, yp, xd, yd, initialSlope, curveConstr)

pp = spline(xp,yp);
yp = ppval(pp,xd); 
F = yp(:)-yd(:);

if(~isempty(initialSlope))
    ppdt = arSplineDer(pp);
    F(end+1) = ppval(ppdt,xp(1)) - initialSlope;
end
if(~isempty(initialSlope))
    if(~exist('ppdt','var'))
        ppdt = arSplineDer(pp);
    end
    ppdtdt = arSplineDer(ppdt);
    tmpF = ppval(ppdtdt,xd) * curveConstr;
    F = [F(:); tmpF(:)];
end