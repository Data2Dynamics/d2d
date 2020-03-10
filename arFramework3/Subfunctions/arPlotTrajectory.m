% [hy, hystd, hyss] = arPlotTrajectory(jy, t, y, ystd, lb, ub, tExp, yExp, yExpHl, yExpStd, ...
%    y_ssa, y_ssa_lb, y_ssa_ub, plotopt, Clines, ClinesExp, trafo, hy, hystd, hyss, dydt, ...
%    qFit[, zero_break, t_ppl, y_ppl_ub, y_ppl_lb, plotOnlyData])
%
% Plots Trajectories

function [hy, hystd, hyss] = arPlotTrajectory(jy, t, y, ystd, lb, ub, tExp, yExp, yExpHl, yExpStd, ...
    y_ssa, y_ssa_lb, y_ssa_ub, plotopt, Clines, ClinesExp, trafo, hy, hystd, hyss, dydt, ...
    qFit, zero_break, t_ppl, y_ppl_ub, y_ppl_lb, plotOnlyData)

fastPlot = false;
if ~exist('plotOnlyData') || isempty(plotOnlyData)
    plotOnlyData = false;
end

% set marker type
if(~isempty(qFit))
    if(qFit(jy))
        ClinesExp{6} = '.';
    else
        ClinesExp{6} = 'o';
    end
end

if ~plotOnlyData
    % plot ssa
    if(~isempty(y_ssa) && any(plotopt(jy)==[3,5]) &&  sum(~isnan(y_ssa))>0)
        y_ssa = trafo(y_ssa);
        
        % create patch around ssa simulation
        if(~isempty(y_ssa_lb))
            y_ssa_lb = trafo(y_ssa_lb);
            y_ssa_ub = trafo(y_ssa_ub);
            for jssa = 1:size(y_ssa_lb, 3)
                tmpx = [t(:); flipud(t(:))];
                tmpy = [y_ssa_ub(:,jy,jssa); flipud(y_ssa_ub(:,jy,jssa))];
                patch(tmpx, tmpy, tmpx*0-2*eps, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2)
            end
        end
        for jssa = 1:size(y_ssa, 3)
            plot(t, y_ssa(:,jy,jssa), 'Color', Clines{2}*0.4+0.6)
        end
        if(size(y_ssa,3)>3)
            plot(t, mean(y_ssa(:,jy,:),3), '--', 'Color', Clines{2})
        end
    end
    
    % Fallback for no data
    if isempty(y)
        hy = plot(NaN, NaN, Clines{:});
        hyss = patch([0 0 0], [0 0 0], ones(1,3));
        set(hyss, 'FaceColor', 'none', 'EdgeColor', 'Black', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2,'LineStyle','none','Marker','*','MarkerEdgeColor',Clines{2},'MarkerSize',5);
        hystd = patch([0 0 0], [0 0 0], -2*ones(1,3), ones(1,3));
        set(hystd, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
        return;
    end
    
    % plot trajectory
    tmpy = y(:,jy);
    tmpy = trafo(tmpy);
    isInfWarn = sum(isinf(tmpy))>0;
    tmpy(isinf(tmpy)) = nan;
    if(isempty(hy))
        hy = plot(t, tmpy, Clines{:});
        %plot data points of model prediction profile likelihood as stars
        if(~isempty(t_ppl) && any(plotopt(jy)==[4,5]))
            hyss = patch([t_ppl(:,jy) ; flipud(t_ppl(:,jy))], [trafo(y_ppl_lb(:,jy)); flipud(trafo(y_ppl_ub(:,jy)))], ones(size([y_ppl_lb(:,jy); y_ppl_ub(:,jy)])));
            set(hyss, 'FaceColor', 'none', 'EdgeColor', 'Black', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2,'LineStyle','none','Marker','*','MarkerEdgeColor',Clines{2},'MarkerSize',5);
        end
        if(any(plotopt(jy) == [3,5]))
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
    if any(plotopt(jy)==[3,4,5])
        tmpx = [t(:); flipud(t(:))];
        tmpy = [];
        if(any(plotopt(jy)==[3,5]) && ~isempty(ystd))
            tmpy = [y(:,jy) + ystd(:,jy); flipud(y(:,jy) - ystd(:,jy))];
        elseif(plotopt(jy)==4 && ~isempty(lb) && ~isempty(ub))
            tmpy = [ub(:,jy); flipud(lb(:,jy))];
        end
        if(~isempty(tmpy))
            tmpy = trafo(tmpy);
            
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
                set(hystd, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2, 'EdgeAlpha', 0.2);
            else
                set(hystd(jy), 'YData', tmpy);
                fastPlot = true;
            end
        end
    end
    
    if(isInfWarn)
        fprintf('Warning: Simulation contains Inf values.\n');
    end
end

% plot data
if(~isempty(qFit) && qFit(jy))
    markersize = 18;
else
    markersize = 6;
end
if(~isempty(yExp) && ~fastPlot)
    yExpU(:,jy)     = trafo(yExp(:,jy)+yExpStd(:,jy)) - trafo(yExp(:,jy));
    yExpL(:,jy)     = trafo(yExp(:,jy)-yExpStd(:,jy)) - trafo(yExp(:,jy));
    yExp            = trafo(yExp);
    yExpHl          = trafo(yExpHl);
    
    if(any(plotopt(jy)==[1,3,4]))
        if(isempty(hy))
            hy = plot(tExp, yExp(:,jy), ClinesExp{:},'MarkerSize',markersize);
        else
            plot(tExp, yExp(:,jy), ClinesExp{:},'MarkerSize',markersize);
        end
        if(sum(~isnan(yExpHl))>0)
            plot(tExp, yExpHl(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',markersize+6);
        end
    elseif any(plotopt(jy)==[2,5])
        errorbar(tExp, yExp(:,jy), yExpL(:,jy), yExpU(:,jy), ClinesExp{:},'MarkerSize',markersize);
        if(sum(~isnan(yExpHl))>0)
            errorbar(tExp, yExpHl(:,jy), yExpL(:,jy), yExpU(:,jy), ClinesExp{:},'LineWidth',2,'MarkerSize',markersize+6);
        end
    end
end

% plot dose response zero line
if(~isempty(zero_break))
    plot([zero_break zero_break], ylim, 'k--');
end


