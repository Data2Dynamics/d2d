% d = arCheckODETolerances(dtol)
% 
%   This function assess the accuracy of ODE intergration by multiplication
%   of atol and rtol with the factors provided as dtol.
%   Default: dtol = 0.1,1,2
% 
%   d.res    maximal difference of ar.res
%   d.sres   maximal difference of ar.sres
%   d.chi2   maximal difference of ar.chi2fit

function d = arCheckODETolerances(dtol)
if(~exist('dtol','var') || isempty(dtol))
    dtol = [.1,1,2];
end

global ar

atolIn = ar.config.atol+0.0;
rtolIn = ar.config.rtol+0.0;

if(~isfield(ar,'sres'))
    arChi2
end
    

sres = NaN(size(ar.sres,1),size(ar.sres,2),length(dtol));
res  = NaN(length(ar.res),length(dtol));
chi2  = NaN(1,length(dtol));
for i=1:length(dtol)
    ar.config.atol = atolIn*dtol(i);
    ar.config.rtol = rtolIn*dtol(i);
    
    try
        arChi2
        res(:,i) = ar.res;
        sres(:,:,i)  = ar.sres;
        chi2(i) = ar.chi2fit;
    catch ERR
        ar.config.atol = atolIn;
        ar.config.rtol = rtolIn;
        rethrow(ERR);
    end
end    

ar.config.atol = atolIn;
ar.config.rtol = rtolIn;

d.res = max(max(abs(diff(res,[],2))));
d.sres = max(max(max(abs(diff(sres,[],3)))));
d.chi2 = max(diff(chi2));


