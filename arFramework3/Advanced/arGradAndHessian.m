% arGradAndHessian
%
% calculate gradients and hessian from sensitivities

function arGradAndHessian

global ar

if ~isfield(ar,'fit')
    warning('arGradAndHessian stopped: ar.fit not available. Fit first.');
    return
end

ar.fit.beta = nan(1, length(ar.p));
ar.fit.alpha = nan(length(ar.p));

ar.fit.beta(ar.fit.qFit==1) = - (ar.fit.res * ar.fit.sres);        % (15.5.6)
ar.fit.gradient = -2 * (ar.fit.beta);               % (15.5.8)

ar.fit.alpha(ar.fit.qFit==1,ar.fit.qFit==1) = ar.fit.sres' * ar.fit.sres;          % (15.5.11)
ar.fit.hessian = 2 * (ar.fit.alpha);                % (15.5.8)

% repress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'all');

ar.fit.covar = nan(length(ar.p));
ar.fit.covar(ar.fit.qFit==1,ar.fit.qFit==1) = inv(ar.fit.alpha(ar.fit.qFit==1,ar.fit.qFit==1));                   % (15.5.15)
% reset warnings
warning(warn_reset);

% eigenvectors and eigenvalues of hessian
ar.fit.hessian_eval = nan(1, length(ar.p));
ar.fit.hessian_evec = nan(length(ar.p));
[evec, eval] = eig(ar.fit.hessian(ar.fit.qFit==1,ar.fit.qFit==1));
ar.fit.hessian_eval(ar.fit.qFit==1) = diag(eval);
ar.fit.hessian_evec(ar.fit.qFit==1,ar.fit.qFit==1) = evec;
