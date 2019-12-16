% arPLEPlot([(jk, singleFigure, do_relative])
% plot sampling of profile likelihood
%
% jk:               parameter index or indices          [all ple calcs]
% singleFigure:     if true, error in code!             [false]
%                                                       [true if jk==1]
% do_relative:      chi2s = chi2s - ar.ple.chi2Reset    [false]

function arPLEPlot(jk, singleFigure, do_relative)

global ar

if(~exist('jk','var') || isempty(jk))
    jk = find(ar.ple.run==1);
end
if(~exist('singleFigure','var') || isempty(singleFigure))
    if(length(jk) == 1)
        singleFigure = true;
    else
        singleFigure = false;
    end
end
if(~exist('do_relative','var'))
    do_relative = false;
end

% iterate over jk in single figures
if(singleFigure)
    if(length(jk)>1)
        for j=1:length(jk)
            arPlotPLE(jk(j));
        end
        return;
    end
end

if(length(jk)==1 && ar.ple.run(jk)==0)
    fprintf('#%i PLE for parameter %s not run yet.\n', jk, ar.pLabel{jk});
    return
end
    
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

if(singleFigure)
    figure(jk)
    nrows = 2;
    ncols = 1;    
else
    figure(1)
    [nrows, ncols] = arNtoColsAndRows(length(jk));
end
clf;

jks = jk;

[~, ~, ylabeltmp] = arGetMerit(true);

gs1 = nan(size(jk));
for j=1:length(jks)
    jk = jks(j);
    chi2s = ar.ple.chi2s{jk};
    if(isfield(ar,'scan'))
        chi2s_scan = ar.scan.chi2s{jk};
    else
        chi2s_scan = [];
    end
    chi2curr = ar.ple.chi2Reset(jk);
    
    if(do_relative)
        chi2curr = chi2curr - ar.ple.chi2Reset(jk);
        chi2s = chi2s - ar.ple.chi2Reset(jk);
    end
    
    ps = ar.ple.ps{jk};
    
    dchi2 = arChi2inv(1-ar.ple.alpha_level, ar.ple.dof_point);
    
    gs1(j) = subplot(nrows, ncols,j);
    
    plot(ps(:,jk), chi2s, 'k');
    arSpacedAxisLimits(gs1(j))
    
    hold on
    if(~isempty(chi2s_scan))
        plot(ps(:,jk), chi2s_scan, 'k--');
    end
    plot(ar.p(jk), chi2curr, '*k');
    
    qerrors = ar.ple.errors{jk} == 0;
    plot(ps(qerrors,jk), chi2s(qerrors), 'ro');
    qerrors = ar.ple.errors{jk} < 0;
    plot(ps(qerrors,jk), chi2s(qerrors), 'rs');
    qerrors = ar.ple.errors{jk} > 1;
%     plot(ps(qerrors,jk), chi2s(qerrors), 'x', 'Color', [.5 .5 .5]);
    plot(ps(qerrors,jk), chi2s(qerrors), 'k.');
    
    % thresholds
    plot(xlim, [0 0]+chi2curr+dchi2, 'r--')
    xlimtmp = xlim;
    text(xlimtmp(1), chi2curr+dchi2, sprintf(' %2i%%', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    
    if(sum(~isnan(chi2s))>0)
        % ylimtmp = [min(chi2s) chi2curr+dchi2*1.2];
        ylimtmp = [min(chi2s) max([max(chi2s)+0.1 chi2curr+dchi2])];
        % ylimtmp = [min([chi2s chi2s_scan]) max([chi2s chi2s_scan])];
        ylim([ylimtmp(1)-diff(ylimtmp)*0.1 ylimtmp(2)+diff(ylimtmp)*0.1]);
    end
    
    % ub and lb
    ylimtmp = ylim;
    if(xlimtmp(1)<ar.lb(jk))
        patch([xlimtmp(1) ar.lb(jk) ar.lb(jk) xlimtmp(1)], [ylimtmp(1) ylimtmp(1) ylimtmp(2) ylimtmp(2) ], ...
            ones(1,4), 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .5);
    end
    if(xlimtmp(2)>ar.ub(jk))
        patch([xlimtmp(2) ar.ub(jk) ar.ub(jk) xlimtmp(2)], [ylimtmp(1) ylimtmp(1) ylimtmp(2) ylimtmp(2)], ...
            ones(1,4), 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .5);
    end
    
    hold off
    
    ylabel(ylabeltmp);
    
    if(singleFigure)
        subplot(2,1,2);
        notjk = (1:length(ar.p));
        notjk = notjk~=jk & ar.qFit==1;
        
        stds = std(ar.ple.ps{jk}(~isnan(ar.ple.chi2s{jk}),notjk));
        [~, istds] = sort(stds, 2, 'descend');
        
        pstmp = ps(:,notjk) - (ones(length(chi2s),1)*ar.ple.pStart(notjk));
        nplot = 7;
        
        if(length(istds)>nplot)
            plot(ps(:,jk), pstmp(:,istds(1:nplot)));
        else
            plot(ps(:,jk), pstmp);
        end
        hold on
        if(length(istds)>nplot)
            plot(ps(:,jk), pstmp(:,istds((nplot+1):end)), 'Color', [.7 .7 .7]);
        end
        plot([ar.ple.pStart(jk) ar.ple.pStart(jk)], ylim, 'k--')
        
        arSpacedAxisLimits(gca)
        
        % ub and lb
        ylimtmp = ylim;
        if(xlimtmp(1)<ar.lb(jk))
            patch([xlimtmp(1) ar.lb(jk) ar.lb(jk) xlimtmp(1)], [ylimtmp(1) ylimtmp(1) ylimtmp(2) ylimtmp(2) ], ...
                ones(1,4), 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .5);
        end
        if(xlimtmp(2)>ar.ub(jk))
            patch([xlimtmp(2) ar.ub(jk) ar.ub(jk) xlimtmp(2)], [ylimtmp(1) ylimtmp(1) ylimtmp(2) ylimtmp(2)], ...
                ones(1,4), 'FaceColor', [.5 .5 .5], 'EdgeColor', 'none', 'FaceAlpha', .5);
        end
        hold off
        
        xlim(xlimtmp);
        ylabel('parameter changes');
        ptmp = ar.pLabel(notjk);
        ptmp = strrep(ptmp,'_','\_');
        if(length(istds)>nplot)
            legend(ptmp{istds(1:nplot)})
        else
            legend(ptmp)
        end
    end
    xlabel(arNameTrafo(ar.pLabel{jk}));
end
% arSpacedAxisLimits(gs1, [], false, true)
