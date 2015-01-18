function [hy, hystd] = arPlotTrajectory(jy, t, y, ystd, lb, ub, tExp, yExp, yExpHl, yExpStd, ...
    y_ssa, y_ssa_lb, y_ssa_ub, ploterrors, Clines, ClinesExp, qUnlog, hy, hystd, qFit, zero_break)

fastPlot = false;

% set marker type
if(qFit)
    ClinesExp{6} = '.';
else
    ClinesExp{6} = 'o';
end

% plot ssa
if(~isempty(y_ssa) && ploterrors==1)
    if(qUnlog)
        y_ssa = 10.^y_ssa;
    end
    
    % create patch around ssa simulation
    if(~isempty(y_ssa_lb))
        if(qUnlog)
            y_ssa_lb = 10.^y_ssa_lb;
            y_ssa_ub = 10.^y_ssa_ub;
        end
        for jssa = 1:size(y_ssa_lb, 3)
            tmpx = [t(:); flipud(t(:))];
            tmpy = [y_ssa_ub(:,jy,jssa); flipud(y_ssa_ub(:,jy,jssa))];
            patch(tmpx, tmpy, tmpx*0-2*eps, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)           
        end
    end
    for jssa = 1:size(y_ssa, 3)
        plot(t, y_ssa(:,jy,jssa), 'Color', Clines{2}*0.4+0.6)
    end
    if(size(y_ssa,3)>3)
        plot(t, mean(y_ssa(:,jy,:),3), '--', 'Color', Clines{2})
    end
end

% plot trajectory
tmpy = y(:,jy);
if(qUnlog)
    tmpy = 10.^tmpy;
end
isInfWarn = sum(isinf(tmpy))>0;
tmpy(isinf(tmpy)) = nan;
if(isempty(hy))
    hy = plot(t, tmpy, Clines{:});
else
    set(hy, 'YData', tmpy);
    fastPlot = true;
end

% plot error bands
if(ploterrors ~= 1)
    tmpx = [t(:); flipud(t(:))];
    tmpy = [];
    if(ploterrors==0 && ~isempty(ystd))
        tmpy = [y(:,jy) + ystd(:,jy); flipud(y(:,jy) - ystd(:,jy))];
    elseif(ploterrors==-1 && ~isempty(lb) && ~isempty(ub))
        tmpy = [ub(:,jy); flipud(lb(:,jy))];
    end
    if(qUnlog)
        tmpy = 10.^tmpy;
    end
    isInfWarn = isInfWarn || sum(isinf(tmpy))>0;
    tmpy(isinf(tmpy)) = nan;
    qnan = isnan(tmpy);
    tmpx = tmpx(~qnan);
    tmpy = tmpy(~qnan);
    
    if(isempty(hystd))
        if(isempty(tmpx))
            hystd = patch([0 0 0], [0 0 0], -2*ones(1,3), ones(1,3));
        else
            hystd = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpy)));
        end
        set(hystd, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    else
        set(hystd, 'YData', tmpy);
        fastPlot = true;
    end
end

if(isInfWarn)
    fprintf('Warning: Simulation contains Inf values.\n');
end

% plot data
if(~isempty(yExp) && ~fastPlot)
    if(qUnlog)
        yExp = 10.^yExp;
        yExpHl = 10.^yExpHl;
    end
    if(ploterrors ~= 1)
        plot(tExp, yExp(:,jy), ClinesExp{:});
        if(sum(~isnan(yExpHl))>0)
            plot(tExp, yExpHl(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',10);
        end
    else
        errorbar(tExp, yExp(:,jy), yExpStd(:,jy), ClinesExp{:});
        if(sum(~isnan(yExpHl))>0)
            errorbar(tExp, yExpHl(:,jy), yExpStd(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',10);
        end
    end
end

% plot dose response zero line
if(~isempty(zero_break))
    plot([zero_break zero_break], ylim, 'k--');
end


