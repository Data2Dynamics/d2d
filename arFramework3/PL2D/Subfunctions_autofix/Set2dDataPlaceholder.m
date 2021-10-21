function [m,d,idpred,iddata] = Set2dDataPlaceholder
% [m,d,idpred,iddata] = Set2dDataPlaceholder
%
% Reads the information on the current 2d-profile and adds a data
% placeholder into the ar-struct which can then be given any value.

global ar

m = ar.ple2d.general.model;
d = ar.ple2d.general.condition;
idpred = ar.ple2d.general.idpred;
arAddToData(m,d,idpred,ar.ple2d.general.tpred,0,2,1);
arLink(true);
iddata = size(ar.model(m).data(d).yExp,1);
ar.model(m).data(d).yExpStd(iddata,idpred) = ar.ple2d.general.sigma;

end

