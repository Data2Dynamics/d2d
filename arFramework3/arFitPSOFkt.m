function f=arFitPSOFkt(p, ar)

try
    ar = arChi2(ar, false,p);
    f = sum(ar.res.^2);
catch
    f = Inf;
end
