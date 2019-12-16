% calculate prediction profile likelihood
%
% arPPL

function arPPL_old(m, c, ix, t, xstd)

global ar

qLog10 = ar.ppl.qLog10;
n = ar.ppl.n;

if(~exist('xstd','var'))
    xstd = 0.01;
end

% optimizer settings (set only once)
if(ar.config.useSensis)
    ar.config.optim.Jacobian = 'on';
else
    ar.config.optim.Jacobian = 'off';
end
ar.config.optim.Algorithm = 'trust-region-reflective';

% chi2 threshold
ar.ppl.dchi2 = arChi2inv(1-ar.ppl.alpha_level, ar.ppl.ndof);

ar.model(m).condition(c).ppl.t = t(:);
ar.model(m).condition(c).ppl.ix = ix(:);
ar.model(m).condition(c).ppl.xtrial = nan(length(t), length(ix), 2*n+1);
ar.model(m).condition(c).ppl.xfit = nan(length(t), length(ix), 2*n+1);
ar.model(m).condition(c).ppl.ppl = nan(length(t), length(ix), 2*n+1);
ar.model(m).condition(c).ppl.ps = nan(length(t), length(ix), 2*n+1, length(ar.p));

ar.model(m).condition(c).ppl.lb_fit = nan(length(t), length(ix));
ar.model(m).condition(c).ppl.ub_fit = nan(length(t), length(ix));

pReset = ar.p + 0;
chi2start = arGetMerit('chi2fit') + 0;    
    
tic;

arWaitbar(0);
ntot = length(t)*length(ix)*n*2;
tcount = 1;
for jt=1:length(t)
    arLink(true, t(jt));
    arChi2(true);
    [~,it] = min(abs(ar.model(m).condition(c).tExp-t(jt)));
    
    for jx = 1:length(ix)
        
        xSim = ar.model(m).condition(c).xExpSimu(it,ix(jx));
        if(qLog10)
            xSim = log10(xSim);
        end
        
        % go up
        npre = (jt-1)*length(ix)*n*2 + (jx-1)*n*2;
        [xtrial_up, xfit_up, ppl_up, ps_up, tcount] = ppl(m, c, ix(jx), it, xSim, ...
            chi2start, 1, qLog10, n, xstd, npre, ntot, tcount);
        
        % go down
        ar.p = pReset;
        arChi2(true);
        npre = (jt-1)*length(ix)*n*2 + (jx-1)*n*2 + n;
        [xtrial_down, xfit_down, ppl_down, ps_down, tcount] = ppl(m, c, ix(jx), it, xSim, ...
            chi2start, -1, qLog10, n, xstd, npre, ntot, tcount);
        
        % reset parameters
        ar.p = pReset;
        arChi2(true);
        
        ar.model(m).condition(c).ppl.xtrial(jt,jx,:) = [fliplr(xtrial_down) xSim xtrial_up];
        ar.model(m).condition(c).ppl.xfit(jt,jx,:) = [fliplr(xfit_down) xSim xfit_up];
        ar.model(m).condition(c).ppl.ppl(jt,jx,:) = [fliplr(ppl_down) chi2start ppl_up];
        ar.model(m).condition(c).ppl.ps(jt,jx,:,:) = [flipud(ps_down); pReset; ps_up];
        
        q_chi2good = squeeze(ar.model(m).condition(c).ppl.ppl(jt,jx,:)) <= chi2start+ar.ppl.dchi2;
        q_nonnan = ~isnan(squeeze(ar.model(m).condition(c).ppl.ppl(jt,jx,:)));
        
        % calculate CI point-wise fit
        ar.model(m).condition(c).ppl.lb_fit(jt,jx) = min(ar.model(m).condition(c).ppl.xfit(jt,jx,q_chi2good));
        ar.model(m).condition(c).ppl.ub_fit(jt,jx) = max(ar.model(m).condition(c).ppl.xfit(jt,jx,q_chi2good));
        if(ar.model(m).condition(c).ppl.lb_fit(jt,jx)==min(ar.model(m).condition(c).ppl.xfit(jt,jx,q_nonnan)))
            ar.model(m).condition(c).ppl.lb_fit(jt,jx) = -Inf;
        else
            kind = find(ar.model(m).condition(c).ppl.xfit(jt,jx,:)==ar.model(m).condition(c).ppl.lb_fit(jt,jx));
            ar.model(m).condition(c).ppl.lb_fit(jt,jx) = interp1(squeeze(ar.model(m).condition(c).ppl.ppl(jt,jx,[kind kind-1])), ...
                squeeze(ar.model(m).condition(c).ppl.xfit(jt,jx,[kind kind-1])), chi2start+ar.ppl.dchi2);
        end
        if(ar.model(m).condition(c).ppl.ub_fit(jt,jx)==max(ar.model(m).condition(c).ppl.xfit(jt,jx,q_nonnan)))
            ar.model(m).condition(c).ppl.ub_fit(jt,jx) = Inf;
        else
            kind = find(ar.model(m).condition(c).ppl.xfit(jt,jx,:)==ar.model(m).condition(c).ppl.ub_fit(jt,jx));
            ar.model(m).condition(c).ppl.ub_fit(jt,jx) = interp1(squeeze(ar.model(m).condition(c).ppl.ppl(jt,jx,[kind kind+1])), ...
                squeeze(ar.model(m).condition(c).ppl.xfit(jt,jx,[kind kind+1])), chi2start+ar.ppl.dchi2);
        end
    end
end
arWaitbar(-1);
toc;

ar.p = pReset;
arChi2(false);

end



function [xtrial, xfit, ppl, ps, tcount] = ppl(m, c, ix, it, xSim, chi2start, direction, qLog10, n, xstd, npre, ntot, tcount)

global ar

xtrial = nan(1,n);
xfit = nan(1,n);
ppl = nan(1,n);
ps = nan(n,length(ar.p));

dx = sqrt(ar.ppl.dchi2*ar.ppl.rel_increase) * xstd;

t = ar.model(m).condition(c).tExp(it);
xLabel = arNameTrafo(ar.model(m).x{ix});

xExp = xSim;
for j = 1:n
    if(toc>tcount)
        if(direction>0)
            arWaitbar((j+npre), ntot, sprintf('PPL (up) for %s at t=%g %i/%i', xLabel, t, j, n));
        else
            arWaitbar((j+npre), ntot, sprintf('PPL (down) for %s at t=%g %i/%i', xLabel, t, j, n));
        end
        tcount = tcount + 0.5; % update every half second
    end

    xExp = xExp + direction*dx;
    
    try
        arPPLFit;
    catch exception
        fprintf('ERROR PPL: going to lower bound (%s)\n', exception.message);
        break;
    end
    
    xtrial(j) = xExp;
    xfit(j) = xSim;
    ppl(j) = arGetMerit('chi2fit');
    ps(j,:) = ar.p;
    
    if(ppl(j) > chi2start+ar.ppl.dchi2*1.2)
        break
    end
end

    function arPPLFit
        
        [pFit, ~, ~, exitflag] = ...
            lsqnonlin(@ppl_merit_fkt, ar.p(ar.qFit==1), ar.lb(ar.qFit==1), ar.ub(ar.qFit==1), ar.config.optim);
        
        ar.p(ar.qFit==1) = pFit;
        
        if(exitflag < 1)
            outputstr = '';
            switch exitflag
                case 1
                    outputstr = 'LSQNONLIN converged to a solution.';
                case 2
                    outputstr = 'Change in X too small.';
                case 3
                    outputstr = 'Change in RESNORM too small.';
                case 4
                    outputstr = 'Computed search direction too small.';
                case 0
                    outputstr = 'Too many function evaluations or iterations.';
                case -1
                    outputstr = 'Stopped by output/plot function.';
                case -2
                    outputstr = 'Bounds are inconsistent.';
                case -3
                    outputstr = 'Regularization parameter too large (Levenberg-Marquardt).';
                case -4
                    outputstr = 'Line search failed.';
            end
            fprintf('lsqnonlin finished after %i interations: %s\n', ar.fit.output.iterations, outputstr);
        end
    end


    function [res, sres] = ppl_merit_fkt(pTrial)
        arChi2(ar.config.useSensis, pTrial)
        
        res = ar.res;
        
        xSim = ar.model(m).condition(c).xExpSimu(it,ix);
        if(qLog10)
            xSim = log10(xSim);
        end
        res(end+1) = (xExp-xSim)/xstd;
        
        if(nargout>1 && ar.config.useSensis)
            sres = ar.sres(:, ar.qFit==1);
            
            sxSim = zeros(1,length(ar.p));
            sxSim(ar.model(m).condition(c).pLink) = ...
                squeeze(ar.model(m).condition(c).sxExpSimu(it,ix,:))' .* ...
                10.^ar.p(ar.model(m).condition(c).pLink) * log(10);
            
            if(qLog10)
                sxSim = sxSim / 10^xSim / log(10);
            end
            
            sres(end+1,:) = -sxSim(ar.qFit==1) / xstd;
        end       
    end

end


function inv = arChi2inv (x, n)
if (nargin ~= 2)
    error ('arChi2inv: you must give two arguments');
end

if (~isscalar (n))
    [retval, x, n] = common_size(x, n);
    if (retval > 0)
        error ('arChi2inv: x and n must be of common size or scalar');
    end
end

inv = gaminv(x, n / 2, 2);
end

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
end
