% Fit only some model parameters ips
%
% arFitSome(ips, silent)

function arFitSome(ips, silent)
global ar

if(nargin<2)
    silent = false;
end

qFitReset = ar.qFit + 0;

qDoFit = ismember(1:length(ar.p),ips);
ar.qFit(ar.qFit==1 & ~qDoFit) = 0;
arFit(silent);

ar.qFit = qFitReset;
