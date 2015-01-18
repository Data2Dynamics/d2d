function [hys, hystds, nrows, ncols] = arPlotTrajectories(ccount, ncount, t, y, ystd, lb, ub, nfine_dr_plot, ...
    nfine_dr_method, tExp, yExp, yExpHl, yExpStd, y_ssa, y_ssa_lb, y_ssa_ub, ...
    ploterrors, qUnlog, qLog, qFit, zero_break, fastPlotTmp, hys, hystds, isLast, qDR, ...
    ndata, chi2, tUnits, response_parameter, yLabel, yNames, yUnits, fiterrors, ...
    logplotting_xaxis)

% rows and cols
ny = size(y,2);
[nrows, ncols] = arNtoColsAndRows(ny);
if(nrows*ncols == ny)
    [nrows, ncols] = arNtoColsAndRows(ny+1);
end

% styles
Clines = arLineMarkersAndColors(ccount, ncount, [], 'none', '-');
ClinesExp = arLineMarkersAndColors(ccount, ncount, [], 'none', 'none');

for jy = 1:ny
    
    % smooth plotted trajectory
    if(length(unique(t))==1)
        t = [t-0.1; t+0.1];
        y = [y; y]; %#ok<AGROW>
        ystd = [ystd; ystd]; %#ok<AGROW>
        lb = [lb; lb]; %#ok<AGROW>
        ub = [ub; ub]; %#ok<AGROW>
    elseif(nfine_dr_plot>10)
        tf = linspace(min(t), max(t), nfine_dr_plot);
        [t, qit] = unique(t);
        y = y(qit);
        ystd = ystd(qit);
        y = interp1(t,y,tf,nfine_dr_method);
        ystd = interp1(t,ystd,tf,nfine_dr_method);
        y = y(:);
        ystd = ystd(:);
        if(~isempty(lb))
            lb = lb(qit);
            ub = ub(qit);
            lb = interp1(t,lb,tf,nfine_dr_method);
            ub = interp1(t,ub,tf,nfine_dr_method);
        end
        t = tf;
    end
    
    if(~fastPlotTmp)
        g = subplot(nrows,ncols,jy);
        hold(g, 'on');
        
        % call arPlotTrajectory
        [hy, hystd] = arPlotTrajectory(jy, t, y, ystd, lb, ub, ...
            tExp, yExp, yExpHl, yExpStd, ...
            y_ssa, y_ssa_lb, y_ssa_ub, ...
            ploterrors, Clines, ClinesExp, qUnlog(jy), ...
            [], [], qFit(jy), zero_break);
        
        % save handles for fast plotting
        hys(jy) = hy;
        if(isempty(hystd))
            hystd = hy;
        end
        hystds(jy) = hystd;
        
        
        % labels and title in last call
        if(isLast)
            arSubplotStyle(g);
            
            qxlabel = jy == (nrows-1)*ncols + 1;
            if(ny <= (nrows-1)*ncols)
                qxlabel = jy == (nrows-2)*ncols + 1;
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
            if(qLog(jy))
                ylabel(g, sprintf('log_{10}(%s) [%s]', yUnits{jy,3}, yUnits{jy,2}));
            else
                ylabel(g, sprintf('%s [%s]', yUnits{jy,3}, yUnits{jy,2}));
            end
            
            titstr = {};
            if(~isempty(yNames) && ~isempty(yNames{jy}) && ~strcmp(yNames{jy}, yLabel{jy}))
                titstr{1} = [arNameTrafo(yNames{jy}) ' (' arNameTrafo(yLabel{jy}) ')'];
            else
                titstr{1} = arNameTrafo(yLabel{jy});
            end
            if(~isempty(yExp))
                if(ndata(jy)>0)
                    if(fiterrors == 1)
                        titstr{2} = sprintf('-2 log(L)_{%i} = %g', ndata(jy), 2*ndata(jy)*log(sqrt(2*pi)) + chi2(jy));
                    else
                        titstr{2} = sprintf('chi^2_{%i} = %g', ndata(jy), chi2(jy));
                    end
                end
            end
            title(g, titstr);
            
            arSpacedAxisLimits(g);
            hold(g, 'off');
        end
    else
        % call arPlotTrajectory
        arPlotTrajectory(jy, t, y, ystd, lb, ub, ...
            tExp, yExp, yExpHl, yExpStd, ...
            y_ssa, y_ssa_lb, y_ssa_ub, ...
            ploterrors, Clines, ClinesExp, qUnlog(jy), ...
            hys(jy), ...
            hystds(jy), qFit(jy), []);
    end
end