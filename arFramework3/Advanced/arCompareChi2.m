function pass = arCompareChi2(ar1,ar2)

if ar1.chi2 == ar2.chi2
    pass = 1;
    return
end
pass = 0;

if ar1.chi2constr ~= ar2.chi2constr
    warning('arCompareChi2.m: ar.chi2constr are different.')
end
if ar1.chi2err ~= ar2.chi2err
    warning('arCompareChi2.m: ar.chi2err are different.')
end
if ar1.chi2fit ~= ar2.chi2fit
    warning('arCompareChi2.m: ar.chi2fit are different.')
end
if ar1.chi2prior ~= ar2.chi2prior
    warning('arCompareChi2.m: ar.chi2prior are different.')
end
if ar1.chi2random ~= ar2.chi2random
    warning('arCompareChi2.m: ar.chi2random are different.')
end