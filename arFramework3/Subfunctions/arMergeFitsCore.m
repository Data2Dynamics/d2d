% ar = arMergeFitsCore(ar,ar2)
% 
% arMergeFitsCore merges the fits in ar2 into ar
% 
%   ar      D2D struct, usually the global ar
%   ar2     Another D2D struct which also contains results of arFits
% 
% 
% See also arFits, arFitLHS, arMergFits

function ar = arMergeFitsCore(ar,ar2)
if ~isstruct(ar) || ~isstruct(ar2)
    error('Both arguments have to be structs')
end
if ~isfield(ar,'ps')
    error('ar does not contain field ar.ps');
end
if ~isfield(ar2,'ps')
    warning('ar2 does not contain field ar2.ps: return without appending.');
end

nfit1 = length(ar.chi2s);
nfitBoth = nfit1 + length(ar2.chi2s);

ar.chi2s((nfit1+1):nfitBoth) = ar2.chi2s;
ar.chi2sconstr((nfit1+1):nfitBoth) = ar2.chi2sconstr;
ar.chi2s_start((nfit1+1):nfitBoth) = ar2.chi2s_start;
ar.chi2sconstr_start((nfit1+1):nfitBoth) = ar2.chi2sconstr_start;
ar.optim_crit((nfit1+1):nfitBoth) = ar2.optim_crit;
ar.ps((nfit1+1):nfitBoth,:) = ar2.ps;
ar.ps_start((nfit1+1):nfitBoth,:) = ar2.ps_start;

ar.exitflag((nfit1+1):nfitBoth) = ar2.exitflag;

