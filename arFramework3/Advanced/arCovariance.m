% arCovariance(silent)
%
% calculate covariance matrix from sensitivities
%
%   silent - boolean if not plotting covariances [true]

function arCovariance(silent)

global ar

if(~exist('silent','var'))
    silent = true;
end

arCalcMerit(true);

ar.alpha = ar.sres' * ar.sres;          % (15.5.11)
ar.hessian = 2 * (ar.alpha);                % (15.5.8)

% repress warning "Matrix is close to singular or badly scaled."
warn_reset = warning;
warning('off', 'all');

ar.covar = nan(size(ar.hessian));
ar.covar(ar.qFit==1,ar.qFit==1) = inv(ar.alpha(ar.qFit==1,ar.qFit==1));                   % (15.5.15)
% reset warnings
warning(warn_reset);

%% plot

if(~silent)
    nroot = 8;
    
    figure(1)
    imagesc(sign(ar.covar).*(abs(ar.covar).^(1/nroot)))
    hold on
    plot(xlim, ylim, 'k:')
    hold off
    colorbar
    set(gca, 'XTickLabel', {});
    set(gca, 'XTick', 1:length(ar.p));
    set(gca, 'YTick', 1:length(ar.p));
    set(gca, 'YTickLabel', ar.pLabel);
    title(sprintf('covariance matrix (scaled by sqrt-%i)',nroot));
end
