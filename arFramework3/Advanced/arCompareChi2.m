function pass = arCompareChi2(ar1,ar2)

fields = {'chi2','chi2constr', 'chi2err', 'chi2fit','chi2prior', 'chi2random'};
pass = 1;
for i = 1:length(fields)
    pass = pass*compareChi2s(ar1, ar2, fields{i});
end
fprintf(['Chi2 of ar1 is ' num2str(ar1.chi2) '. Chi2 of ar2 is ' num2str(ar2.chi2) '.\n'])

end

function pass = compareChi2s(ar1, ar2, field)

%atol = 1e-3;
rtol = 1e-2;

absComp = abs(ar1.(field) - ar2.(field));
if abs(ar1.(field)) ~= 0
    relComp = absComp/abs(ar1.(field));
elseif abs(ar1.(field)) + abs(ar2.(field)) == 0
    relComp = 0;
else 
    relComp = Inf;
end

if ar1.(field) == ar2.(field)
    pass = 1;
elseif  relComp < rtol % && absComp < atol   % eg Borghans: ar1.chi2 = 8445.5 ar2.chi2 = 8445.6 due to log/lin data 
    fprintf(['arCompareChi2: ar.' field ' are different within tolerances.\n'])
    pass = 1;
else
    warning(['arCompareChi2.m: ar.' field ' are different.'])
    pass = 0;
end

end