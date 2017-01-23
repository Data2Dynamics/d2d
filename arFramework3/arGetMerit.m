%   [meritval, meritvals. meritLabel] = arGetMerit
%   [meritval, meritvals. meritLabel] = arGetMerit(silent)
%   [meritval, meritvals. meritLabel] = arGetMerit(whichone)
% 
%   silent      false/true, default: false
%   whichone    character indicating which merit should be calculated
%               default: '' (then the optimization merit is calculated)
%               'lsqnonlin'
%               'loglik_all': the total likelihood incl. ar.chi2constr
%               'loglik': the likelihood except ar.chi2constr
%               'chi2fit': total chi2 except chi2constr
%               'chi2_all': the total chi2 incl. chi2constr
% 
%               The argument whichone can be used to keep the old
%               terminology.
% 
% 
%   This function prints different merits/objective functions which can
%   be used for optimization/parameter estimation.
% 
%   ATTENTION: arCalcMerit, arSimu, ... are NOT called automatically (only
%   if ar.res_type is missing). To update residuals and merits this has to
%   be done independently. 
% 
%       meritval    1) ~isempty(whichone): the requested merit
% 
%                   2) meritval = '': value of the merit function currently
%                   used for parameter fitting 
% 
%       meritvals   different optional merits and sub-terms
% 
%       meritLabel  name/label of meritval as it can be used for plot
%       labels
% 
%   See: 
%   https://github.com/Data2Dynamics/d2d/wiki/Objective%20function,%20likelhood%20and%20chi-square%20in%20the%20d2d%20framework
% 
% 
%   Example 1: Printing merit at the command line:
%   arGetMerit
% 
%   Example 2: Value of objective function currently used for optimization: 
%   objfun = arGetMerit
% 
%   Example 3: Value of objective function currently used for optimization
%   and other merits and merit terms:
%   [meritval, allmerits] = arGetMerit
% 
%   Example 4: silent:
%   meritval = arGetMerit(true);
% 
%   Example 5: Special merit value:
%   chi2prior = arGetMerit('chi2prior');

function varargout = arGetMerit(arg1)
if nargin==1
    if ischar(arg1)
        whichone = arg1;
        silent = true;  % if only a special merit is requested, then always silent
    elseif islogical(arg1) || isnumeric(arg1)
        silent = arg1;        
    end
end

if ~exist('silent','var') || isempty(silent)
    silent = false;
end
if ~exist('whichone','var') || isempty(whichone)
    whichone = '';  % empty means all merits
end


global ar 

if isempty(ar)
    error('global ar is empty, please load a model.');
elseif ~isfield(ar,'res_type')
    arCalcMerit(true,ar.p(ar.qFit==1)); %sensi = true ensures that ODE intergation steps are the same as within fitting
end


%% 1) Fast calls if nargout==1 
% fast access to chi2fit, chi2, ar.chi2err, chi2prior and chi2constr 
if nargout == 1
    switch whichone
        case 'chi2fit'
            varargout{1} = ar.chi2fit;
            return
        case 'chi2'
            varargout{1} = ar.chi2;
            return
        case 'chi2err'
            varargout{1} = ar.chi2err;
            return
        case 'chi2prior'
            varargout{1} = ar.chi2prior;
            return
        case 'chi2constr'
            varargout{1} = ar.chi2constr;
            return
    end
end

%% 2) Other function calls
meritvals.loglik        = 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit;
meritvals.loglik_all    = 2*(ar.ndata+ar.nconstr)*log(sqrt(2*pi)) + ar.chi2fit + ar.chi2constr;
meritvals.lsqnonlin     = sum([ar.res,ar.constr].^2);

meritvals.chi2_all      = ar.chi2fit + ar.chi2constr;
meritvals.chi2_res      = sum(ar.res(ar.res_type==1).^2); % same as ar.chi2 - ar.chi2prior
meritvals.chi2_err      = ar.chi2err;
meritvals.chi2_err_addc = sum(ar.res(ar.res_type==2).^2); % sum of error residuals
meritvals.chi2_constr   = ar.chi2constr;
meritvals.chi2_prior    = ar.chi2prior;
meritvals.chi2_random   = ar.chi2random;

if ~isempty(ar.ple) && isfield(ar.ple,'merit_fkt')    
    meritvals.ple_merit = feval(ar.ple.merit_fkt);
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


%% meritval
switch lower(whichone)  % case insensitive
    case 'lsqnonlin'
        meritval = meritvals.lsqnonlin;
    
    case 'loglik'
        meritval = meritvals.loglik;
        meritLabel = '-2 log likelihood'
    
    case {'loglik_all','ll','likelihood','loglik_tot'}
        meritval = meritvals.loglik_all;
        meritLabel = '-2 log likelihood total';
        
    case {'chi2all','chi2_all'}
        meritval = meritvals.chi2_all;
        meritLabel = '\chi^2_{total}';

    case 'chi2fit'
        meritval = ar.chi2fit;
        meritLabel = '\chi^2_{fit}';

    case 'chi2prior'
        meritval = meritvals.chi2_prior;
        meritLabel = '\chi^2_{prior}';

    case 'chi2random'
        meritval = meritvals.chi2_random;
        meritLabel = '\chi^2_{random}';

    case 'chi2constr'
        meritval = meritvals.chi2_constr;
        meritLabel = '\chi^2_{constr}';

    case 'chi2'
        meritval = meritvals.chi2;
        meritLabel = '\chi^2';
        
    otherwise  % default merit used for optimization
        if ~isempty(whichone)
            error('arGetMerit.m: whichone=%s is not implemented',whichone)
        else
            if ar.config.fiterrors == -1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)==0) % if no error parameters fitted
                meritval = meritvals.chi2_all;
                meritLabel = '\chi^2_{total}';
                if(~silent)
                    arFprintf(1, 'total chi^2 = %g, ', meritvals.chi2_all);
                end
            else
                meritval = meritvals.loglik_all;
                meritLabel = 'log likelihood';
                if(~silent)
                    arFprintf(1, '-2*log(L) = %g, %i data points, ', ...
                        meritvals.loglik, meritvals.ndata);
                end
            end
        end
end

%% output arguments
if nargout>0
    varargout{1} = meritval;
end

if nargout>1
    varargout{2} = meritvals;
end

if nargout>2
    varargout{3} = meritLabel;
end


%% printing at the command line
if(~silent)
    arFprintf(1, '%i free parameters, ', sum(ar.qFit==1));
    if(ar.chi2constr ~=0)
        arFprintf(1, ', %g violation of %i constraints', ar.chi2constr, ar.nconstr);
    end
    if(ar.chi2prior ~=0)
        arFprintf(1, ', %g violation of %i prior assumptions', ar.chi2prior, ar.nprior);
    end
    arFprintf(1, 'data chi^2 = %g', meritvals.chi2_res);
%     if(sensi)
%         arFprintf(1, ', first order optimality criterion %g (%i)', ar.firstorderopt, -sum(qred));
%     end
    arFprintf(1, '\n');
end


    