function [hy, hystd] = arPlotTrajectory(g, t, y, ystd, lb, ub, tExp, yExp, yExpHl, yExpStd, ...
    y_ssa, y_ssa_lb, y_ssa_ub, ploterrors, Clines, ClinesExp, qUnlog)

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
            tmpy = [y_ssa_ub(:,1,jssa); flipud(y_ssa_ub(:,1,jssa))];
            patch(tmpx, tmpy, tmpx*0-2*eps, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
            hold(g, 'on');
            
        end
    end
    for jssa = 1:size(y_ssa, 3)
        plot(t, y_ssa(:,1,jssa), 'Color', Clines{2}*0.4+0.6)
        hold(g, 'on');
    end
    if(size(y_ssa,3)>3)
        plot(t, mean(y_ssa(:,1,:),3), '--', 'Color', Clines{2})
    end
end

tmpy = y;
if(qUnlog)
    tmpy = 10.^tmpy;
end
isInfWarn = sum(isinf(tmpy))>0;
tmpy(isinf(tmpy)) = nan;
hy = plot(g, t, tmpy, Clines{:});

hold(g, 'on');
if(ploterrors ~= 1)
    tmpx = [t(:); flipud(t(:))];
    tmpy = [];
    if(ploterrors==0 && ~isempty(ystd))
        tmpy = [y + ystd; flipud(y - ystd)];
    elseif(ploterrors==-1 && ~isempty(lb) && ~isempty(ub))
        tmpy = [ub; flipud(lb)];
    end
    if(qUnlog)
        tmpy = 10.^tmpy;
    end
    isInfWarn = isInfWarn || sum(isinf(tmpy))>0;
    tmpy(isinf(tmpy)) = nan;
    qnan = isnan(tmpy);
    tmpx = tmpx(~qnan);
    tmpy = tmpy(~qnan);
    
    if(isempty(tmpx))
        ltmp = patch([0 0 0], [0 0 0], -2*ones(1,3), ones(1,3));
    else
        ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpy)));
    end
    set(ltmp, 'FaceColor', Clines{2}, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    hystd = ltmp;
end

if(isInfWarn)
    fprintf('Warning: Simulation contains Inf values.\n');
end

if(~isempty(yExp))
    if(qUnlog)
        yExp = 10.^yExp;
        yExpHl = 10.^yExpHl;
    end
    if(ploterrors ~= 1)
        plot(g, tExp, yExp, ClinesExp{:});
        if(sum(~isnan(yExpHl))>0)
            hold(g,'on');
            plot(g, tExp, yExpHl, ClinesExp{:},'LineWidth',2,'MarkerSize',10);
        end
    else
        errorbar(g, tExp, yExp, yExpStd, ClinesExp{:});
        if(sum(~isnan(yExpHl))>0)
            hold(g,'on');
            errorbar(g, tExp, yExpHl, yExpStd, ClinesExp{:},'LineWidth',2,'MarkerSize',10);
        end
    end
end