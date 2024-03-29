% Plot 2D log-Likelihood for arFitNonLin
%
% plePlot2D(jk1, jk2, [figure_number, nsteps, nlevels])
%
% figure_number     [1]
% nsteps            [100]
% nlevels           [40]

function plePlot2D(jk1, jk2, figure_number, nsteps, nlevels)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end

if(~ar.qFit(jk1) || ~ar.qFit(jk2))
    error('specified parameters need to be free to fit');
end

if(~exist('figure_number', 'var'))
    figure_number = 1;
end
if(~exist('nsteps', 'var'))
    nsteps = 100;
end
if(~exist('nlevels', 'var'))
    nlevels = 40;
end

%% sampling range
oversamplefactor = 0.1;

k1low = min(ar.ple.conf_lb(:,jk1));
k1high = max(ar.ple.conf_ub(:,jk1));
if(isnan(k1low) || k1low==-Inf)
    k1low = min(ar.ple.ps{jk1}(:,jk1));
end
if(isnan(k1high) || k1high==Inf)
    k1high = max(ar.ple.ps{jk1}(:,jk1));
end
k1range = k1high-k1low;
k1steps = linspace(k1low-k1range*oversamplefactor, k1high+k1range*oversamplefactor, nsteps);

k2low = min(ar.ple.conf_lb(:,jk2));
k2high = max(ar.ple.conf_ub(:,jk2));
if(isnan(k2low) || k2low==-Inf)
    k2low = min(ar.ple.ps{jk2}(:,jk2));
end
if(isnan(k2high) || k2high==Inf)
    k2high = max(ar.ple.ps{jk2}(:,jk2));
end
k2range = k2high-k2low;
k2steps = linspace(k2low-k2range*oversamplefactor, k2high+k2range*oversamplefactor, nsteps);

%% sample
chi2change = nan(length(k1steps),length(k2steps));

preset = ar.ple.p;
h = waitbar(0, 'Please wait...');
for j1=1:length(k1steps)
    waitbar(j1/length(k1steps), h);
    for j2=1:length(k2steps)
        ptmp = preset;
        ptmp([jk1 jk2]) = [k1steps(j1) k2steps(j2)];
        feval(ar.ple.integrate_fkt, ptmp);
        chi2change(j2, j1) = feval(ar.ple.merit_fkt);
    end
end
feval(ar.ple.integrate_fkt, preset);
close(h)

%% plot
chi2levels_amplification = 10;

dchi2change1D = ar.ple.dchi2_point;
dchi2change = ar.ple.dchi2;

figure(figure_number)

% plot contours
chi2levels = logspace(log10(min(chi2change(:))*1.01), ...
    log10(min(chi2change(:)) + dchi2change*chi2levels_amplification), nlevels);
contour(k1steps,k2steps,chi2change, chi2levels);
% mycmap = sqrt(gray(500));
mycmap = gray(500);
colormap(mycmap)
colorbar
hold on

% plot CI
contour(k1steps,k2steps,chi2change,min(chi2change(:))+dchi2change1D, ...
    'Color', [.5 0 0], 'LineWidth', 2);
contour(k1steps,k2steps,chi2change,min(chi2change(:))+dchi2change, ...
    'Color', [1 0 0], 'LineWidth', 2);

% plot truth
if(isfield(ar.ple, 'p_true'))
    plot(ar.ple.p_true(jk1), ar.ple.p_true(jk2), 'g*')
end

% plot PLE
if(isfield(ar.ple, 'ps'))
    if(length(ar.ple.ps)>=jk1)
        tmp_ps = ar.ple.ps{jk1};
        tmp_psstart = ar.ple.psinit{jk1};
        if(~isempty(tmp_ps))
            tmp_chi2s = ar.ple.chi2s{jk1};
            q_tmp_chi2s = tmp_chi2s < min(tmp_chi2s)+ar.ple.dchi2;
            q_tmp_chi2s_point = tmp_chi2s < min(tmp_chi2s)+ar.ple.dchi2_point;
            plot(tmp_ps(:,jk1), tmp_ps(:,jk2), 'b--')
            plot(tmp_ps(q_tmp_chi2s,jk1), tmp_ps(q_tmp_chi2s,jk2), 'b-')
            plot(tmp_ps(q_tmp_chi2s_point,jk1), tmp_ps(q_tmp_chi2s_point,jk2), 'b-', ...
                'LineWidth', 2)
            plot(tmp_psstart(:,jk1), tmp_psstart(:,jk2), 'rx')
        end
    end
    if(length(ar.ple.ps)>=jk2)
        tmp_ps = ar.ple.ps{jk2};
        tmp_psstart = ar.ple.psinit{jk2};
        if(~isempty(tmp_ps))
            tmp_chi2s = ar.ple.chi2s{jk2};
            q_tmp_chi2s = tmp_chi2s < min(tmp_chi2s)+ar.ple.dchi2;
            q_tmp_chi2s_point = tmp_chi2s < min(tmp_chi2s)+ar.ple.dchi2_point;
            plot(tmp_ps(:,jk1), tmp_ps(:,jk2), 'b--')
            plot(tmp_ps(q_tmp_chi2s,jk1), tmp_ps(q_tmp_chi2s,jk2), 'b-')
            plot(tmp_ps(q_tmp_chi2s_point,jk1), tmp_ps(q_tmp_chi2s_point,jk2), 'b-', ...
                'LineWidth', 2)
            plot(tmp_psstart(:,jk1), tmp_psstart(:,jk2), 'rx')
        end
    end
end

% plot hessian
if(sum(isnan(ar.ple.conf_lb(1,jk1)) & isnan(ar.ple.conf_ub(1,jk1))) == 0)
    plot(ar.ple.p(jk1), ar.ple.p(jk2), '*r')
    try
        error_ellipse(ar.ple.covar([jk1 jk2], [jk1 jk2]), [ar.ple.p(jk1) ar.ple.p(jk2)], 'conf', 1-ar.ple.alpha_level)
    catch
        fprintf('Covariance Matrix not positive-semi-definit!!!\n');
    end
    % plot([0 0] + ar.ple.p(jk1) - ar.ple.conf_lb(1,jk1), ylim, 'b--')
    % plot([0 0] + ar.ple.p(jk1) + ar.ple.conf_ub(1,jk1), ylim, 'b--')
    % plot(xlim, [0 0] + ar.ple.p(jk2) - ar.ple.conf_lb(1,jk2), 'b--')
    % plot(xlim, [0 0] + ar.ple.p(jk2) + ar.ple.conf_ub(1,jk2), 'b--')
end


hold off
title(sprintf('\\chi^2 and %5.2f%% confidence regions', (1-ar.ple.alpha_level)*100))
if(ar.qLog10(jk1))
    xlabel(['log_{10}(' ar.ple.p_labels{jk1} ')'])
else
    xlabel(ar.ple.p_labels{jk1})
end
if(ar.qLog10(jk2))
    ylabel(['log_{10}(' ar.ple.p_labels{jk2} ')'])
else
    ylabel(ar.ple.p_labels{jk2})
end

% save
if(exist(ar.ple.savePath, 'dir'))
    saveas(gcf, [ar.ple.savePath '/2D_' ar.ple.p_labels{jk1} '_' ar.ple.p_labels{jk2}], 'fig')
    saveas(gcf, [ar.ple.savePath '/2D_' ar.ple.p_labels{jk1} '_' ar.ple.p_labels{jk2}], 'eps')
end

