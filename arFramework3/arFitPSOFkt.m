function f=arFitPSOFkt(p, ar)

global fit

fit.fevals = fit.fevals + 1;
try
    ar = arChi2(ar, false,p);
    f = ar.chi2fit;
catch
    f = Inf;
end
