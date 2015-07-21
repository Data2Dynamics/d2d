function [hy, hystd, hyss] = arPlotTrajectory(jy, t, y, ystd, lb, ub, tExp, yExp, yExpHl, yExpStd, ...
    y_ssa, y_ssa_lb, y_ssa_ub, ploterrors, Clines, ClinesExp, qUnlog, qLog, hy, hystd, hyss, dydt, ...
    qFit, zero_break, t_ppl, y_ppl_ub, y_ppl_lb)

fastPlot = false;

% set marker type
if(~isempty(qFit))
    if(qFit(jy))
        ClinesExp{6} = '.';
    else
        ClinesExp{6} = 'o';
    end
end

% plot ssa
if(~isempty(y_ssa) && ploterrors==1 &&  sum(~isnan(y_ssa))>0)
    if(~isempty(qUnlog) && qUnlog(jy))
        y_ssa = 10.^y_ssa;
    end
    if(~isempty(qLog) && qLog(jy))
        y_ssa = log10(y_ssa);
    end
    
    % create patch around ssa simulation
    if(~isempty(y_ssa_lb))
        if(~isempty(qUnlog) && qUnlog(jy))
            y_ssa_lb = 10.^y_ssa_lb;
            y_ssa_ub = 10.^y_ssa_ub;
        end
        if(~isempty(qLog) && qLog(jy))
            y_ssa_lb = log10(y_ssa_lb);
            y_ssa_ub = log10(y_ssa_ub);
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
if(~isempty(qUnlog) && qUnlog(jy))
    tmpy = 10.^tmpy;
end
if(~isempty(qLog) && qLog(jy))
    tmpy = log10(tmpy);
end
isInfWarn = sum(isinf(tmpy))>0;
tmpy(isinf(tmpy)) = nan;
if(isempty(hy))
    hy = plot(t, tmpy, Clines{:});
    %plot data points of model prediction profile likelihood as stars
    if(~isempty(t_ppl))
        if(~isempty(qUnlog) && qUnlog(jy))
            hyss = patch([t_ppl(:,jy) ; flipud(t_ppl(:,jy))], [10.^y_ppl_lb(:,jy); flipud(10.^y_ppl_ub(:,jy))], ones(size([y_ppl_lb(:,jy); y_ppl_ub(:,jy)])));             
        else
            hyss = patch([t_ppl(:,jy) ; flipud(t_ppl(:,jy))], [y_ppl_lb(:,jy); flipud(y_ppl_ub(:,jy))], ones(size([y_ppl_lb(:,jy); y_ppl_ub(:,jy)])));             
        end
       set(hyss, 'FaceColor', 'none', 'EdgeColor', 'Black', 'FaceAlpha', 0.2);
    end
    if(ploterrors ~= 1)
        set(hy, 'LineWidth', 1.5);
    end 
else
    set(hy(jy), 'YData', tmpy);
    fastPlot = true;
end

% steady state
if(~isempty(dydt) && length(dydt)>=jy)
    if(~isnan(dydt(jy)))
        yss = y(1,jy) + dydt(jy)*(t-min(t));
        if(isempty(hyss))
            ClinesSS = Clines;
            ClinesSS{4} = '--';
            hyss = plot(t, yss, ClinesSS{:});
        else
            set(hyss(jy), 'YData', yss);
        end
    else
        hyss = nan;
    end
else
    hyss = [];
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
    if(~isempty(tmpy))
        if(~isempty(qUnlog) && qUnlog(jy))
            tmpy = 10.^tmpy;
        end
        if(~isempty(qLog) && qLog(jy))
            tmpy = log10(tmpy);
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
            set(hystd(jy), 'YData', tmpy);
            fastPlot = true;
        end
    end
end

if(isInfWarn)
    fprintf('Warning: Simulation contains Inf values.\n');
end

% plot data
if(~isempty(qFit) && qFit(jy))
    markersize = 18;
else
    markersize = 6;
end
if(~isempty(yExp) && ~fastPlot)
    if(~isempty(qUnlog) && qUnlog(jy))
        yExp = 10.^yExp;
        yExpHl = 10.^yExpHl;
    end
    if(~isempty(qLog) && qLog(jy))
        yExp = log10(yExp);
        yExpHl = log10(yExpHl);
    end
    if(ploterrors ~= 1)
        plot(tExp, yExp(:,jy), ClinesExp{:},'MarkerSize',markersize);
        if(sum(~isnan(yExpHl))>0)
            plot(tExp, yExpHl(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',markersize+6);
        end
    else
        errorbar(tExp, yExp(:,jy), yExpStd(:,jy), ClinesExp{:},'MarkerSize',markersize);
        if(sum(~isnan(yExpHl))>0)
            errorbar(tExp, yExpHl(:,jy), yExpStd(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',markersize+6);
        end
    end
end

% plot dose response zero line
if(~isempty(zero_break))
    plot([zero_break zero_break], ylim, 'k--');
end


