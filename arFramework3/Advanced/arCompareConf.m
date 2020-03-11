function pass = arCompareConf(ar1,ar2)

rtol = 1e-3;
pass = 1;

if ~(ar1.config.rtol-ar2.config.rtol==0)
    pass = 0;
    warning('arCompareConf.m: ar.config.rtol is different.')
end
if ~(ar1.config.atol-ar2.config.atol==0)
    pass = 0;
    warning('arCompareConf.m: ar.config.atol is different.')
end
if max( ( (ar1.res-ar2.res) ./ ar1.res ) .^2 )>rtol
    pass = 0;
    warning('arCompareConf.m: ar.res are different.')
end
if ar1.config.fiterrors~=ar2.config.fiterrors
    warning('arCompareConf.m: ar.config.fiterrors are set different.');
end