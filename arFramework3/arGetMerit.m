%   [meritval, meritvals] = arGetMerit
% 
%   This function prints different merits/objective functions which can
%   be used for optimization/parameter estimation.
% 
%   ATTENTION: arCalcMerit, arSimu, ... are NOT called automatically (only
%   if ar.res_type is missing). To update residuals and merits this has to
%   be done independently. 
% 
%       meritval    value of the merit function currently used for parameter fitting
% 
%       meritvals   different optional merits and sub-terms
% 
%   See: 
%   https://github.com/Data2Dynamics/d2d/wiki/Objective%20function,%20likelhood%20and%20chi-square%20in%20the%20d2d%20framework
% 
% 
%   Example 1: Printing merit at the command line:
% 
%   arGetMerit
% 
% 
%   Example 2: Value of objective function currently used for optimization:
% 
%   objfun = arGetMerit
% 
% 
%   Example 3: Value of objective function currently used for optimization
%   and other merits and merit terms:
% 
%   [objfun, allmerits] = arGetMerit

function varargout = arGetMerit(silent)

if ~exist('silent','var') || isempty(silent)
    silent = false;
end


global ar pleGlobals

if ~isfield(ar,'res_type')
    arCalcMerit(true,ar.p(ar.qFit==1)); %sensi = true ensures that ODE intergation steps are the same as within fitting
end

meritvals.loglik = 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit;
meritvals.lsqnonlin = sum([ar.res,ar.constr].^2);

meritvals.chi2_all      = ar.chi2fit;
meritvals.chi2_res      = sum(ar.res(ar.res_type==1).^2); % same as ar.chi2 - ar.chi2prior
meritvals.chi2_err      = ar.chi2err;
meritvals.chi2_err_addc = sum(ar.res(ar.res_type==2).^2);
meritvals.chi2_constr   = ar.chi2constr;
meritvals.chi2_prior    = ar.chi2prior;
meritvals.chi2_random   = ar.chi2random;

if ~isempty(pleGlobals)
    meritvals.ple_merit = feval(pleGlobals.merit_fkt);
else
    meritvals.ple_merit = [];
end

meritvals.nres  = length(ar.res(:));
meritvals.ndata = ar.ndata;
meritvals.nerr = sum(ar.res_type==2);
meritvals.nconstr = ar.nconstr;
meritvals.nprior  = ar.nprior;
meritvals.nconstr = ar.nconstr;

meritvals.npara = sum(ar.qFit==1);

meritvals.fiterrors = ar.config.fiterrors;
meritvals.fiterrors_correction = ar.config.fiterrors_correction;


if ar.config.fiterrors == -1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)==1)==0) % if no error parameters fitted
    meritval = meritvals.chi2_res;
else
    meritval = meritvals.loglik;
end

if nargout>0
    varargout{1} = meritval;
end

if nargout>1
    varargout{2 } = meritvals;
end

if(~silent)
    arFprintf(1, '-2*log(L) = %g, %i data points, ', ...
        meritvals.loglik, meritvals.ndata);
    arFprintf(1, '%i free parameters,  ', sum(ar.qFit==1));
    arFprintf(1, 'total chi^2 = %g, ', meritvals.chi2_all);
    arFprintf(1, 'data chi^2 = %g', meritvals.chi2_res);
    if(ar.chi2constr ~=0)
        arFprintf(1, ', %g violation of %i constraints', ar.chi2constr, ar.nconstr);
    end
    if(ar.chi2prior ~=0)
        arFprintf(1, ', %g violation of %i prior assumptions', ar.chi2prior, ar.nprior);
    end
%     if(sensi)
%         arFprintf(1, ', first order optimality criterion %g (%i)', ar.firstorderopt, -sum(qred));
%     end
    arFprintf(1, '\n');
end


    