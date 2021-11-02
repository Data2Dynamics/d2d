function chi2 = UsePLEFit(idpar)

% chi2 = UsePLEFit(idpar)
%
% Returns objective function value for objective function which is used for
% normal PLE-function


global ar

ar.qFit(idpar) = 0;
arFit(true);
ar.qFit(idpar) = 1;
arCalcMerit(ar.ple.usesensis, ar.p(ar.qFit==1));
chi2 = arGetMerit(true);

end

