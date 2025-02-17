% d = arCheckODETolerances(dtol)
% 
%   This function assess the accuracy of ODE intergration by multiplication
%   of atol and rtol with the factors provided as dtol.
%   Default: dtol = 0.1,1,2
% 
%   d.res    maximal difference of ar.res
%   d.sres   maximal difference of ar.sres
%   d.chi2   maximal difference of arGetMerit('chi2fit')
% 
% Example:
% arCheckODETolerances
% 
% Example:
% arCheckODETolerances(logspace(-2,2,5))

function [d,res,sres] = arCheckODETolerances(dtol)
if(~exist('dtol','var') || isempty(dtol))
    dtol = [.1,1,2];
end

global ar 

atolIn = ar.config.atol+0.0;
rtolIn = ar.config.rtol+0.0;

if(~isfield(ar,'sres'))
    arCalcMerit
end
    

sres = NaN(size(ar.sres,1),size(ar.sres,2),length(dtol));
res  = NaN(length(ar.res),length(dtol));
chi2  = NaN(1,length(dtol));
pleMerit  = NaN(1,length(dtol));
for i=1:length(dtol)
    ar.config.atol = atolIn*dtol(i);
        ar.config.rtol = rtolIn*dtol(i);

    fprintf('atol=%.3e, rtol=%.3e\n',ar.config.atol,ar.config.rtol);
    try
        arCalcMerit(true,ar.p(ar.qFit==1))
        res(:,i) = ar.res;
        sres(:,:,i)  = ar.sres;
        chi2(i) = arGetMerit('chi2fit');
        if ~isempty(ar.ple) && isfield(ar.ple,'merit_fkt')
            pleMerit(i) = feval(ar.ple.merit_fkt);
        end
    catch ERR
        ar.config.atol = atolIn;
        ar.config.rtol = rtolIn;
        if ~isempty(regexp(ERR.message,'cvodes failed at CV_TOO_MUCH_WORK', 'once'))
            warning(ERR.message);
        elseif ~isempty(regexp(ERR.message,'NaN in derivative of residuals', 'once'))            
            warning(ERR.message);
        elseif ~isempty(regexp(ERR.message,'CV_CONV_FAILURE', 'once'))            
            warning(ERR.message);
        else
            rethrow(ERR);
        end
    end
end    

ar.config.atol = atolIn;
ar.config.rtol = rtolIn;

d.res = max(max(abs(range(res,2))));
d.sres = max(max(max(abs(range(sres,3)))));
d.chi2 = max(range(chi2));
d.pleMerit = max(range(pleMerit));

fprintf('Maximal absolute impact of the tolerances:\n');
fprintf('%20s\t%25s:\t %e\n','chi2','(arGetMerit(''chi2fit''))',d.chi2);
fprintf('%20s\t%25s:\t %e\n','pleMerit','(ar.ple.merit_fkt)',d.pleMerit);
fprintf('%20s\t%25s:\t %e\n','Residuals','(ar.res)',d.res);
fprintf('%20s\t%25s:\t %e\n','Sensititivites','(ar.sres)',d.sres);

