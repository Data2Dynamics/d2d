% Initialize Profile Likelihood Exploit
%
% pleInit(p, q_fit, lb, ub, hb_dist, q_log10, integrate_fkt, objective_fkt, diff_fkt, ...
%         fit_fkt, setoptim_fkt, covar_fkt[, p_labels, stdalpha])
%
% p:                         parameter vector
% q_fit:                     which parameters are fitted
% lb:                        lower parameter bound
% ub:                        upper parameter bound
% q_log10:                   which parameters are logarithmic (basis 10)
% integrate_fkt(p):          set parameters and update objective function
% chi2 = objective_fkt:      get merit value
% [beta, alpha] = diff_fkt:  calculate beta=-0.5*dchi^2/dp and alpha=0.5*d^2chi^2/dp^2
% [p, covar] = fit_fkt(i):   call optimization with i'th parameter fixed
% setoptim_fkt(p):           set persistent optimal parameter values
% p_labels:                  cell array of parameter labels
% stdalpha:                  confidence level in standard deviation [1] = 68%

function pleInit(p, q_fit, lb, ub, q_log10, integrate_fkt, merit_fkt, ...
    diffmerit_fkt, fit_fkt, setoptim_fkt, p_labels, alpha, force)
if(~exist('force','var') || isempty(force))
    force = true;
end

global pleGlobals;
if(~isstruct(pleGlobals) || force)
    pleGlobals = struct('p', p, 'q_fit', q_fit, 'lb', lb, 'ub', ub, ...
        'q_log10', q_log10, ...
        'integrate_fkt', integrate_fkt, 'merit_fkt', merit_fkt, ...
        'diffmerit_fkt', diffmerit_fkt, 'fit_fkt', fit_fkt, ...
        'setoptim_fkt', setoptim_fkt, ...
        'initstep_fkt', @pleInitStepDirect, ...  % is overwritten in arPLEInit, if mode is provided
        'dof', sum(q_fit), 'dof_point', 1);
end

pleGlobals(1).allowbetteroptimum = false;

% step sizes
pleGlobals.samplesize = 50 * ones(size(pleGlobals.p));
pleGlobals.relchi2stepincrease = 0.1 * ones(size(pleGlobals.p));
pleGlobals.maxstepsize = (ub-lb)./pleGlobals.samplesize;
pleGlobals.minstepsize = ones(size(pleGlobals.maxstepsize))*1e-3;
pleGlobals.breakonlb = false(size(pleGlobals.p));
pleGlobals.breakonub = false(size(pleGlobals.p));
pleGlobals.attempts = 4;

pleGlobals.showCalculation = true;

pleGlobals.ylabel = '\chi^2_{PL}';

pleGlobals.breakon_point = true; %false;
pleGlobals.plot_point = true;
pleGlobals.plot_simu = false; %true;

pleGlobals.dist_thres = 0.01;
pleGlobals.grad_thres = 1;

pleGlobals.closetobound = 0.001;

% alpha level
if(~exist('alpha', 'var'))
    alpha = 0.95;
end
pleGlobals.alpha_level = 1-alpha;

% delta chi^2 thresholds
pleGlobals.dchi2 = chi2inv(1-pleGlobals.alpha_level, pleGlobals.dof);
pleGlobals.dchi2_point = chi2inv(1-pleGlobals.alpha_level, pleGlobals.dof_point);

% magic factors
pleGlobals.chi2_strID_ratio = 1e-1;
pleGlobals.svd_threshold = 1e-6; % SVD regulatization threshold (NR: chapter 15.4)
pleGlobals.optimset_tol = 1e-1;

% labels
if(~exist('p_labels', 'var'))
    p_labels = {};
    for j=1:length(p)
        p_labels{j} = sprintf('p%02i', j); %#ok<AGROW>
    end
end
pleGlobals.p_labels = p_labels;

pleGlobals.conf_lb = nan(1,length(p));
pleGlobals.conf_ub = nan(1,length(p));
pleGlobals.conf_lb_point = nan(1,length(p));
pleGlobals.conf_ub_point = nan(1,length(p));
pleGlobals.conf_rel = nan(1,length(p));
pleGlobals.conf_rel_point = nan(1,length(p));

pleGlobals.IDstatus = nan(1,length(p));
pleGlobals.IDstatus_point = nan(1,length(p));
pleGlobals.IDlabel = {'', 'pra.nID', 'str.nID', 'single str.nID'};
pleGlobals.savePath = ['PLE-' datestr(now, 30)];

if(~isfield(pleGlobals,'ps') || force)
    pleGlobals.ps = {};
end
if(~isfield(pleGlobals,'psinit') || force)
    pleGlobals.psinit = {};
end
if(~isfield(pleGlobals,'psinitstep') || force)
    pleGlobals.psinitstep = {};
end
if(~isfield(pleGlobals,'chi2s') || force)
    pleGlobals.chi2s = {};
end
if(~isfield(pleGlobals,'chi2sinit') || force)
    pleGlobals.chi2sinit = {};
end

pleGlobals.finished = 0;

% fprintf('\nProfile Likelihood Approach for Identifiability Analysis\n\n');
% fprintf('For detailed information please read:\n');
% fprintf('A. Raue, C. Kreutz, T. Maiwald, J. Bachmann, M. Schilling, U. Klingmueller and J. Timmer\n');
% fprintf('Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood.\n');
% 
% fprintf('Bioinformatics (2009), 25(15), 1923-1929.\n');
% site = '"http://bioinformatics.oxfordjournals.org/cgi/content/abstract/btp358v1"';
% title = 'doi:10.1093/bioinformatics/btp358';
% fprintf(['<a href = ' site '>' title '</a>'])
% fprintf('\n\n');

function inv = chi2inv (x, n)
if (nargin ~= 2)
    error ('chi2inv: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size(x, n);
    if (retval > 0)
        error ('chi2inv: x and n must be of common size or scalar');
    end
end

inv = gaminv(x, n / 2, 2);

function inv = gaminv (x, a, b)
if (nargin ~= 3)
    error ('gaminv: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gaminv: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
inv = zeros (sz);

k = find ((x < 0) | (x > 1) | isnan (x) | ~(a > 0) | ~(b > 0));
if (any (k))
    inv (k) = NaN;
end

k = find ((x == 1) & (a > 0) & (b > 0));
if (any (k))
    inv (k) = Inf;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (~isscalar(a) || ~isscalar(b))
        a = a (k);
        b = b (k);
        y = a .* b;
    else
        y = a * b * ones (size (k));
    end
    x = x (k);
    l = find (x < eps);
    if (any (l))
        y(l) = sqrt (eps) * ones (length (l), 1);
    end
    
    y_old = y;
    for i = 1 : 100

        h     = (gamcdf (y_old, a, b) - x) ./ gampdf (y_old, a, b);
        y_new = y_old - h;
        ind   = find (y_new <= eps);
        if (any (ind))
            y_new (ind) = y_old (ind) / 10;
            h = y_old - y_new;
        end
        if (max (abs (h)) < sqrt (eps))
            break;
        end
        y_old = y_new;
    end
    
    inv (k) = y_new;
end

function pdf = gampdf (x, a, b)
if (nargin ~= 3)
    error ('gampdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gampdf: x, a and b must be of common size or scalars');
    end
end

sz = size(x);
pdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    pdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (a <= 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = (x(k) .^ (a - 1)) ...
            .* exp(- x(k) ./ b) ./ gamma (a) ./ (b .^ a);
    else
        pdf(k) = (x(k) .^ (a(k) - 1)) ...
            .* exp(- x(k) ./ b(k)) ./ gamma (a(k)) ./ (b(k) .^ a(k));
    end
end

k = find ((x > 0) & (a > 1) & (b > 0));
if (any (k))
    if (isscalar(a) && isscalar(b))
        pdf(k) = exp (- a .* log (b) + (a-1) .* log (x(k)) ...
            - x(k) ./ b - gammaln (a));
    else
        pdf(k) = exp (- a(k) .* log (b(k)) + (a(k)-1) .* log (x(k)) ...
            - x(k) ./ b(k) - gammaln (a(k)));
    end
end




function cdf = gamcdf (x, a, b)
if (nargin ~= 3)
    error ('gamcdf: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('gamcdf: x, a and b must be of common size or scalars');
    end
end

sz = size (x);
cdf = zeros (sz);

k = find (~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    cdf (k) = NaN;
end

k = find ((x > 0) & (a > 0) & (b > 0));
if (any (k))
    if (isscalar (a) && isscalar(b))
        cdf (k) = gammainc (x(k) ./ b, a);
    else
        cdf (k) = gammainc (x(k) ./ b(k), a(k));
    end
end
