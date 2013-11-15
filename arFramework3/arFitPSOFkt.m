function f=arFitPSOFkt(p, ~)

global ar

arChi2(false,p);
f = sum(ar.res.^2);
