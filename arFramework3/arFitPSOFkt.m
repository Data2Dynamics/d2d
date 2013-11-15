function f=arFitPSOFkt(p, ~)

global ar

try
    arChi2(false,p);
    f = sum(ar.res.^2);
catch
    f = Inf;
end
