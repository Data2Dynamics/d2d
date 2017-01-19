function [hys, hystds, hysss, nrows, ncols] = arPlotTrajectories(ccount, ncount, t, y, ystd, lb, ub, nfine_dr_plot, ...
    nfine_dr_method, tExp, yExp, yExpHl, yExpStd, y_ssa, y_ssa_lb, y_ssa_ub, ...
    plotopt, qUnlog, qLog, qLogPlot, qFit, zero_break, fastPlotTmp, hys, hystds, hysss, dydt, isLast, qDR, ...
    ndata, chi2, tUnits, response_parameter, yLabel, yNames, yUnits, fiterrors, ...
    logplotting_xaxis, iy, t_ppl, y_ppl_ub, y_ppl_lb)

% rows and cols
ny = length(iy);
[nrows, ncols] = arNtoColsAndRows(ny);
if(nrows*ncols==ny && nrows*ncols>0)
    [nrows, ncols] = arNtoColsAndRows(ny+1);
end 

% styles
Clines = arLineMarkersAndColors(ccount, ncount, [], 'none', '-');
ClinesExp = arLineMarkersAndColors(ccount, ncount, [], 'none', 'none');

for jys = 1:length(iy)
    jy = iy(jys);
    
    % smooth plotted trajectory
    if(length(unique(t))==1)
        t = [t-0.1; t+0.1];
        y = [y; y]; %#ok<AGROW>
        ystd = [ystd; ystd]; %#ok<AGROW>
        lb = [lb; lb]; %#ok<AGROW>
        ub = [ub; ub]; %#ok<AGROW>
    elseif(nfine_dr_plot>10 && nfine_dr_plot>2*length(t))
        tf = linspace(min(t), max(t), nfine_dr_plot)';
        [t, qit] = unique(t);
        y = y(qit,:);
        y = interp1(t,y,tf,nfine_dr_method);
        if(~isempty(ystd))
            ystd = ystd(qit,:);
            ystd = interp1(t,ystd,tf,nfine_dr_method);
        end
        if(~isempty(lb))
            lb = lb(qit,:);
            ub = ub(qit,:);
            lb = interp1(t,lb,tf,nfine_dr_method);
            ub = interp1(t,ub,tf,nfine_dr_method);
        end
        t = tf;
    end
    
    if(~fastPlotTmp)
        g = subplot(nrows,ncols,jys);
        hold(g, 'on');
        
        % call arPlotTrajectory
        [hy, hystd, hyss] = arPlotTrajectory(jy, t, y, ystd, lb, ub, ...
            tExp, yExp, yExpHl, yExpStd, ...
            y_ssa, y_ssa_lb, y_ssa_ub, ...
            plotopt, Clines, ClinesExp, qUnlog, qLog, ...
            [], [], [], dydt, qFit, zero_break, t_ppl, y_ppl_ub, y_ppl_lb);
        
        % save handles for fast plotting
        hys(jy) = hy;
        if(isempty(hystd))
            hystd = hy;
        end
        hystds(jy) = hystd;
        if(~isempty(hyss))
            hysss(jy) = hyss;
        end
        
        % labels and title in last call
        if(isLast)
            arSubplotStyle(g);
            
            % labels
            qxlabel = jys == (nrows-1)*ncols + 1;
            if(ny <= (nrows-1)*ncols)
                qxlabel = jys == (nrows-2)*ncols + 1;
            end
            if(qxlabel)
                if(~qDR)
                    xlabel(g, sprintf('%s [%s]', tUnits{3}, tUnits{2}));
                else
                    if(logplotting_xaxis)
                        xlabel(g, sprintf('log_{10}(%s)', arNameTrafo(response_parameter)));
                    else
                        xlabel(g, sprintf('%s', arNameTrafo(response_parameter)));
                    end
                end
            end
            yunitetmp = '';
            if(~isempty(yUnits{jy,2}))
                yunitetmp = sprintf(' [%s]', yUnits{jy,2});
            end
            if(qLogPlot(jy))
                ylabel(g, sprintf('log_{10}(%s)%s', yUnits{jy,3}, yunitetmp));
            else
                ylabel(g, sprintf('%s%s', yUnits{jy,3}, yunitetmp));
            end
            
            % title
            if(~isempty(yNames) && ~isempty(yNames{jy}) && ~strcmp(yNames{jy}, yLabel{jy}))
                titstr = [arNameTrafo(yNames{jy}) ' (' arNameTrafo(yLabel{jy}) ')'];
            else
                titstr = arNameTrafo(yLabel{jy});
            end
            title(g, titstr);
            
            % text
            if(~isempty(yExp))
                if(ndata(jy)>0)
                    if fiterrors == 1
                        titstr = sprintf('  -2 log(L)_{%i} = %g', ndata(jy), 2*ndata(jy)*log(sqrt(2*pi)) + chi2(jy));
                    else
                        titstr = sprintf('  chi^2_{%i} = %g', ndata(jy), chi2(jy));
                    end
                    text(.95,.95, titstr, 'Units', 'normalized', 'FontSize', 8, ...
                        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
                end
            end
            
            arSpacedAxisLimits(g);
            hold(g, 'off');
        end
    else
        % call arPlotTrajectory
        arPlotTrajectory(jy, t, y, ystd, lb, ub, ...
            tExp, yExp, yExpHl, yExpStd, ...
            y_ssa, y_ssa_lb, y_ssa_ub, ...
            plotopt, Clines, ClinesExp, qUnlog, qLog, ...
            hys, hystds, hysss, dydt, qFit, []);
        arSpacedAxisLimits(get(hys(jy),'Parent'));
    end
end
