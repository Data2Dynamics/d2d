% Plot models Y
%
% arPlotY(saveToFile, fastPlot)
%
% saveToFile    [false]
% fastPlot      [false]
% doLegends      [true]

function arPlotY(saveToFile, fastPlot, doLegends)

matVer = ver('MATLAB');

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('saveToFile','var'))
	saveToFile = false;
end
if(~exist('fastPlot','var'))
	fastPlot = false;
end
if(~exist('doLegends','var'))
	doLegends = true;
end

% constants
labelfontsize = 12;
labelfonttype = 'Arial';
rowstocols = 0.5;
overplot = 0.1;

logplotting_xaxis = true;

figcount = 1;
for jm = 1:length(ar.model)
    ar.model(jm).chi2 = 0;
    ar.model(jm).ndata = 0;
    
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotYs(jplot)==1 && ar.model(jm).plot(jplot).ny>0)
            if(ar.config.ploterrors == -1)
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['CI-Y: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            else
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['Y: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            end
            
            % plotting
            ccount = 1;
            chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
            ndata = zeros(1,ar.model(jm).plot(jplot).ny);   
            if(~ar.model(jm).plot(jplot).doseresponse)
                cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                
                for jd = ar.model(jm).plot(jplot).dLink
                    % rows and cols
                    [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols);
                    
                    Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount);
                    
                    for jy = 1:ny
                        [t, y, ystd, tExp, yExp, yExpStd, lb, ub] = getData(jm, jd, jy);
                        if(~fastPlotTmp)
                            g = subplot(nrows,ncols,jy);
                            ar.model(jm).plot(jplot).gy(jy) = g;
                            
                            if(ar.model(jm).data(jd).qFit(jy))
                                markerstyle = '*';
                            else
                                markerstyle = 'o';
                            end
                            
                            if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                % plot ssa
                                if(isfield(ar.model(jm).data(jd), 'yFineSSA') && ar.config.ploterrors==1)
                                    if(isfield(ar.model(jm).data(jd), 'yFineSSA_lb'))
                                        for jssa = 1:size(ar.model(jm).data(jd).yFineSSA_lb, 3)
                                            tmpx = [t(:); flipud(t(:))];
                                            tmpy = [10.^ar.model(jm).data(jd).yFineSSA_ub(:,jy,jssa); ...
                                                flipud(10.^ar.model(jm).data(jd).yFineSSA_lb(:,jy,jssa))];
                                            patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
%                                              patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.4+0.6, 'FaceColor', Clines{2}*0.4+0.6)
                                            hold(g, 'on');
                                        end
                                    end
                                    for jssa = 1:size(ar.model(jm).data(jd).yFineSSA, 3)
                                        plot(t, 10.^ar.model(jm).data(jd).yFineSSA(:,jy,jssa), 'Color', Clines{2}*0.4+0.6)
                                        hold(g, 'on');
                                    end
                                    if(size(ar.model(jm).data(jd).yFineSSA,3)>1)
                                        plot(t, 10.^mean(ar.model(jm).data(jd).yFineSSA(:,jy,:),3), '--', 'Color', Clines{2})
                                    end
                                end
                                
                                ar.model(jm).data(jd).plot.y(jy) = plot(g, t, 10.^y, Clines{:});
                                cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                hold(g, 'on');
                                if(ar.config.ploterrors ~= 1)
                                    tmpx = [t(:); flipud(t(:))];
                                    if(ar.config.ploterrors==0)
                                        tmpy = [10.^(y + ystd); flipud(10.^(y - ystd))];
                                    elseif(ar.config.ploterrors==-1)
                                        tmpy = [10.^ub; flipud(10.^lb)];
                                    end
                                    ltmp = patch(tmpx, tmpy, -2*ones(size(tmpx)), ones(size(tmpx)));
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, -ones(size(tmpx)), ones(size(tmpx)));
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                    ar.model(jm).data(jd).plot.ystd(jy) = ltmp;
                                    ar.model(jm).data(jd).plot.ystd2(jy) = ltmp2;
                                end
                                
                                if(isfield(ar.model(jm).data(jd), 'yExp'))
                                    if(ar.config.ploterrors ~= 1)
                                        plot(g, tExp, 10.^yExp, markerstyle, Clines{:});
                                    else
                                        errorbar(g, tExp, 10.^yExp, ...
                                            10.^yExp - 10.^(yExp - yExpStd), 10.^(yExp + yExpStd) - 10.^yExp, markerstyle, Clines{:});
                                    end
                                end
                            else
                                % plot ssa
                                if(isfield(ar.model(jm).data(jd), 'yFineSSA') && ar.config.ploterrors==1)
                                    if(isfield(ar.model(jm).data(jd), 'yFineSSA_lb'))
                                        for jssa = 1:size(ar.model(jm).data(jd).yFineSSA_lb, 3)
                                            tmpx = [t(:); flipud(t(:))];
                                            tmpy = [ar.model(jm).data(jd).yFineSSA_ub(:,jy,jssa); ...
                                                flipud(ar.model(jm).data(jd).yFineSSA_lb(:,jy,jssa))];
                                            patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
                                            hold(g, 'on');
                                            
                                        end
                                    end
                                    for jssa = 1:size(ar.model(jm).data(jd).yFineSSA, 3)
                                        plot(t, ar.model(jm).data(jd).yFineSSA(:,jy,jssa), 'Color', Clines{2}*0.4+0.6)
                                        hold(g, 'on');
                                    end
                                    if(size(ar.model(jm).data(jd).yFineSSA,3)>1)
                                        plot(t, mean(ar.model(jm).data(jd).yFineSSA(:,jy,:),3), '--', 'Color', Clines{2})
                                    end
                                end
                                
                                tmpx = t;
                                tmpy = y;
                                qfinite = ~isinf(tmpy);
                                if(sum(qfinite)>0)
                                    ar.model(jm).data(jd).plot.y(jy) = plot(g, tmpx(qfinite), tmpy(qfinite), Clines{:});
                                    cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                end
                                hold(g, 'on');
                                if(ar.config.ploterrors ~= 1)
                                    tmpx = [t(:); flipud(t(:))];
                                    if(ar.config.ploterrors==0)
                                        tmpy = [y + ystd; flipud(y - ystd)];
                                    elseif(ar.config.ploterrors==-1)
                                        tmpy = [ub; flipud(lb)];
                                    end
                                    qfinite = ~isinf(tmpy);
                                    if(sum(qfinite)>0)
                                        ltmp = patch(tmpx(qfinite), tmpy(qfinite), -2*ones(size(tmpx(qfinite))), ones(size(tmpy(qfinite))));
                                        set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                        ltmp2 = patch(tmpx(qfinite), tmpy(qfinite), -ones(size(tmpx(qfinite))), ones(size(tmpy(qfinite))));
                                        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                        ar.model(jm).data(jd).plot.ystd(jy) = ltmp;
                                        ar.model(jm).data(jd).plot.ystd2(jy) = ltmp2;
                                    end
                                end
                                
                                if(isfield(ar.model(jm).data(jd), 'yExp'))
                                    if(ar.config.ploterrors ~= 1)
                                        plot(g, tExp, yExp, markerstyle, Clines{:});
                                    else
                                        errorbar(g, tExp, yExp, yExpStd, markerstyle, Clines{:});
                                    end
                                end
                            end
                        else
                            if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                set(ar.model(jm).data(jd).plot.y(jy), 'YData', 10.^y);
                                if(ar.config.fiterrors ~= -1 && ar.config.ploterrors ~= 1)
                                    set(ar.model(jm).data(jd).plot.ystd(jy), 'YData', [10.^(y + ystd); flipud(10.^(y - ystd))]);
                                    set(ar.model(jm).data(jd).plot.ystd2(jy), 'YData', [10.^(y + ystd); flipud(10.^(y - ystd))]);
                                end
                            else
                                tmpy = y;
                                qfinite = ~isinf(tmpy);
                                set(ar.model(jm).data(jd).plot.y(jy), 'YData', tmpy(qfinite));
                                if(ar.config.fiterrors ~= -1 && ar.config.ploterrors ~= 1)
                                    tmpy = [y + ystd; flipud(y - ystd)];
                                    qfinite = ~isinf(tmpy);
                                    if(sum(qfinite)>0)
                                        set(ar.model(jm).data(jd).plot.ystd(jy), 'YData', tmpy(qfinite));
                                        set(ar.model(jm).data(jd).plot.ystd2(jy), 'YData', tmpy(qfinite));
                                    end
                                end
                            end
                        end
                        
                        % chi^2 & ndata
                        if(ar.model(jm).data(jd).qFit(jy)==1 && ~isempty(ar.model(jm).data(jd).chi2))
                            chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                            ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                            if(ar.config.fiterrors==1)
                                chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2err(jy);
                            end
                        end
                    end
                    ccount = ccount + 1;
                end
            else
                times = [];
                for jd = ar.model(jm).plot(jplot).dLink
					times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                    [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols);

                    for jy = 1:ny
                        % chi^2 & ndata
                        if(ar.model(jm).data(jd).qFit(jy)==1)
                            chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                            ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                            if(ar.config.fiterrors==1)
                                chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2err(jy);
                            end
                        end
                    end
                end
                
                if(str2double(matVer.Version)>=8.1)
                    [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition,'legacy'); %#ok<ASGLU>
                else
                    [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition); %#ok<ASGLU>
                end
                
                cclegendstyles = zeros(1,length(times)*length(conditions));
                for jt = 1:length(times)
                    if(isempty(conditions))
                        jcs = 1;
                    else
                        jcs = 1:length(conditions);
                    end
                    for jc = jcs
                        if(isempty(conditions))
                            ds = ar.model(jm).plot(jplot).dLink;
                        else
                            ds = ar.model(jm).plot(jplot).dLink(find(jconditions==jc)); %#ok<FNDSB>
                        end
                        
                        jd = ds(1);
                        Clines = myLineStyle(length(times)*length(jcs), ccount);
                        
                        for jy = 1:ny
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break] = ...
                                getDataDoseResponse(jm, jy, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                            
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,jy);
                                ar.model(jm).plot(jplot).gy(jy) = g;
                                
                                markerstyle = '*';
                                
                                if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                    qfinite = ~isinf(t) & ~isinf(y);
                                    
                                    ar.model(jm).data(jd).plot.y(jy,jt,jc) = plot(g, t(qfinite), 10.^y(qfinite), Clines{:});
                                    cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy,jt,jc);
                                    hold(g, 'on');
                                    if(ar.config.ploterrors ~= 1)
                                        tmpx = [t(:); flipud(t(:))];
                                        if(ar.config.ploterrors==0)
                                            tmpy = [10.^(y + ystd); flipud(10.^(y - ystd))];
                                        elseif(ar.config.ploterrors==-1)
                                            tmpy = [10.^ub; flipud(10.^lb)];
                                        end
                                        qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                        ltmp = patch(tmpx(qfinite), tmpy(qfinite), tmpx(qfinite)*0-2*eps, 'r');
                                        set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                        ltmp2 = patch(tmpx(qfinite), tmpy(qfinite), tmpx(qfinite)*0-eps, 'r');
                                        set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                        ar.model(jm).data(jd).plot.ystd(jy,jt,jc) = ltmp;
                                        ar.model(jm).data(jd).plot.ystd2(jy,jt,jc) = ltmp2;
                                    end
                                    if(isfield(ar.model(jm).data(jd), 'yExp'))
                                        if(ar.config.ploterrors~=1)
                                            plot(g, tExp, 10.^yExp, markerstyle, Clines{:});
                                        else
                                            errorbar(g, tExp, 10.^yExp, 10.^yExp - 10.^(yExp - yExpStd), ...
                                                10.^(yExp + yExpStd) - 10.^yExp, markerstyle, Clines{:});
                                        end
                                    end
                                else
                                    tmpx = t;
                                    tmpy = y;
                                    qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                    ar.model(jm).data(jd).plot.y(jy,jt,jc) = plot(g, tmpx(qfinite), tmpy(qfinite), Clines{:});
                                    cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy,jt,jc);
                                    hold(g, 'on');
                                    if(ar.config.ploterrors ~= 1)
                                        tmpx = [t(:); flipud(t(:))];
										if(ar.config.ploterrors==0)
											tmpy = [y + ystd; flipud(y - ystd)];
										elseif(ar.config.ploterrors==-1)
											tmpy = [ub; flipud(lb)];
										end
                                        qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                        if(sum(qfinite)>0)
                                            ltmp = patch(tmpx(qfinite), tmpy(qfinite), -2*eps*ones(size(tmpy(qfinite))), 'r');
                                            set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                            ltmp2 = patch(tmpx(qfinite), tmpy(qfinite), -eps*ones(size(tmpy(qfinite))), 'r');
                                            set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                            ar.model(jm).data(jd).plot.ystd(jy,jt,jc) = ltmp;
                                            ar.model(jm).data(jd).plot.ystd2(jy,jt,jc) = ltmp2;
                                        end
                                    end
                                    if(isfield(ar.model(jm).data(jd), 'yExp'))
                                        if(ar.config.ploterrors~=1)
                                            plot(g, tExp, yExp, markerstyle, Clines{:});
                                        else
                                            errorbar(g, tExp, yExp, yExpStd, markerstyle, Clines{:});
                                        end
                                    end
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                    qfinite = ~isinf(t) & ~isinf(y);
                                    set(ar.model(jm).data(jd).plot.y(jy,jt,jc), 'YData', 10.^y(qfinite));
                                    if(ar.config.fiterrors ~= -1)
                                        tmpx = [t(:); flipud(t(:))];
                                        if(ar.config.ploterrors==0)
                                            tmpy = [10.^(y + ystd); flipud(10.^(y - ystd))];
                                        elseif(ar.config.ploterrors==-1)
                                            tmpy = [10.^ub; flipud(10.^lb)];
                                        end
                                        qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                        set(ar.model(jm).data(jd).plot.ystd(jy,jt,jc), 'YData', tmpy(qfinite));
                                        set(ar.model(jm).data(jd).plot.ystd2(jy,jt,jc), 'YData', tmpy(qfinite));
                                    end
                                else
                                    tmpx = t;
                                    tmpy = y;
                                    qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                    set(ar.model(jm).data(jd).plot.y(jy,jt,jc), 'YData', tmpy(qfinite));
                                    if(ar.config.fiterrors ~= -1)
                                        tmpx = [t(:); flipud(t(:))];
										if(ar.config.ploterrors==0)
											tmpy = [y + ystd; flipud(y - ystd)];
										elseif(ar.config.ploterrors==-1)
											tmpy = [ub; flipud(lb)];
										end
                                        qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                        if(sum(qfinite)>0)
                                            set(ar.model(jm).data(jd).plot.ystd(jy,jt,jc), 'YData', tmpy(qfinite));
                                            set(ar.model(jm).data(jd).plot.ystd2(jy,jt,jc), 'YData', tmpy(qfinite));
                                        end
                                    end
                                end
                            end
                        end
                        ccount = ccount + 1;
                    end
                end
            end
            
            % axis & titles
            for jy = 1:ny
                g = ar.model(jm).plot(jplot).gy(jy);
                if(~fastPlotTmp)
                    hold(g, 'off');
                    mySubplotStyle(g, labelfontsize, labelfonttype);
                    
                    if(jy == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).tUnits{3}, ar.model(jm).data(jd).tUnits{2}));
                        else
                            if(logplotting_xaxis)
                                xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                            else
                                xlabel(g, sprintf('%s', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                            end
                        end
                    end
                    if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
                        ylabel(g, sprintf('log_{10}(%s) [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    else
                        ylabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    end
                    
                    if(doLegends && jy == 1 && (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            if(length(ar.model(jm).plot(jplot).dLink)>1)
                                legend(g, cclegendstyles, myNameTrafo(ar.model(jm).plot(jplot).condition))
                            end
                        else
                            legendtmp = {};
                            ccount = 1;
                            for jt=1:length(times)
                                if(~isempty(conditions))
                                    for jc = 1:length(conditions)
                                        legendtmp{ccount} = sprintf('t=%g%s : %s', times(jt), ar.model(jm).tUnits{2}, conditions{jc}); %#ok<AGROW>
                                        ccount = ccount + 1;
                                    end
                                else
                                    legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                    ccount = ccount + 1;
                                end
                            end
                            legend(g, cclegendstyles, myNameTrafo(legendtmp))
                        end
                    end
                end
                titstr = {};
                if(isfield(ar.model(jm).data(jd), 'yNames') && ~isempty(ar.model(jm).data(jd).yNames{jy}) && ...
                        ~strcmp(ar.model(jm).data(jd).yNames{jy}, ar.model(jm).data(jd).y{jy}))
                    titstr{1} = [myNameTrafo(ar.model(jm).data(jd).yNames{jy}) ' (' myNameTrafo(ar.model(jm).data(jd).y{jy}) ')']; %#ok<AGROW>
                else
                    titstr{1} = myNameTrafo(ar.model(jm).data(jd).y{jy}); %#ok<AGROW>
                end
                if(isfield(ar.model(jm).data(jd), 'yExp'))
                    if(ndata(jy)>0)
                        if(ar.config.fiterrors == 1)
                            titstr{2} = sprintf('-2 log(L)_{%i} = %g', ndata(jy), 2*ndata(jy)*log(sqrt(2*pi)) + chi2(jy)); %#ok<AGROW>
                        else
                            titstr{2} = sprintf('chi^2_{%i} = %g', ndata(jy), chi2(jy)); %#ok<AGROW>
                        end
                    end
                end
                title(g, titstr);
                arSpacedAxisLimits(g, overplot);
            end
            
            ar.model(jm).plot(jplot).chi2 = sum(chi2);
            ar.model(jm).plot(jplot).ndata = sum(ndata);
            
            ar.model(jm).chi2 = ar.model(jm).chi2 + sum(chi2);
            ar.model(jm).ndata = ar.model(jm).ndata + sum(ndata);
            
            figcount = figcount + 1;
            
            if(saveToFile)
                if(ar.config.ploterrors == -1)
                    ar.model(jm).plot(jplot).savePath_FigYCI = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                else
                    ar.model(jm).plot(jplot).savePath_FigY = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                end
            end
        else
            try %#ok<TRYNC>
                close(ar.model(jm).plot(jplot).fighandel_y)
            end
            ar.model(jm).plot(jplot).fighandel_y = [];
        end
    end
    
    
end



function [t, y, ystd, tExp, yExp, yExpStd, lb, ub] = getData(jm, jd, jy)
global ar

t = ar.model(jm).data(jd).tFine;
y = ar.model(jm).data(jd).yFineSimu(:,jy);
ystd = ar.model(jm).data(jd).ystdFineSimu(:,jy);
if(isfield(ar.model(jm).data(jd), 'yExp') && ~isempty(ar.model(jm).data(jd).yExp))
    tExp = ar.model(jm).data(jd).tExp;
    yExp = ar.model(jm).data(jd).yExp(:,jy);
    if(ar.config.fiterrors == -1)
        yExpStd = ar.model(jm).data(jd).yExpStd(:,jy);
    else
        yExpStd = ar.model(jm).data(jd).ystdExpSimu(:,jy);
    end
else
    tExp = [];
    yExp = [];
	yExpStd = [];
end
if(isfield(ar.model(jm).data(jd), 'yFineLB'))
	lb = ar.model(jm).data(jd).yFineLB(:,jy);
	ub = ar.model(jm).data(jd).yFineUB(:,jy);
else
	lb = [];
	ub = [];
end



function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break] = getDataDoseResponse(jm, jy, ds, ttime, dLink, logplotting_xaxis)
global ar

zero_break = [];

ccount = 1;
for jd = ds
	qt = ar.model(jm).data(jd).tExp == ttime;
    for jt = find(qt')
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(1).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(1).value); %#ok<AGROW>
        end
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(1).value))))
                        doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(1).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(1).value); %#ok<AGROW>
                end
            end
			doses = unique(doses); %R2013a compatible
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        tExp(ccount,1) = t(ccount,1); %#ok<AGROW>
        
        [tdiffmin, jtfine] = min(abs(ar.model(jm).data(jd).tFine-ar.model(jm).data(jd).tExp(jt))); %#ok<ASGLU>
        y(ccount,1) = ar.model(jm).data(jd).yFineSimu(jtfine,jy); %#ok<AGROW>
        ystd(ccount,1) = ar.model(jm).data(jd).ystdFineSimu(jtfine,jy); %#ok<AGROW>
        
        yExp(ccount,1) = ar.model(jm).data(jd).yExp(jt,jy); %#ok<AGROW>
        if(ar.config.fiterrors == -1)
            yExpStd(ccount,1) = ar.model(jm).data(jd).yExpStd(jt,jy); %#ok<AGROW>
        else
            yExpStd(ccount,1) = ar.model(jm).data(jd).ystdExpSimu(jt,jy); %#ok<AGROW>
        end
        if(isfield(ar.model(jm).data(jd), 'yExpUB'))
            lb(ccount,1) = ar.model(jm).data(jd).yFineLB(jtfine,jy); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).data(jd).yFineUB(jtfine,jy); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end



%% sub-functions

function C = myLineStyle(n, j)
farben = lines(n);
farben(1,:) = [0 0 0];
C = cell(1,2);
C{1} = 'Color';
C{2} = farben(j,:);



function [h, fastPlotTmp] = myRaiseFigure(m, jplot, figname, figcount, fastPlot)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;
fastPlotTmp = fastPlot;

if(ar.config.ploterrors == -1)
    if(isfield(ar.model(m).plot(jplot), 'fighandel_yCI') && ~isempty(ar.model(m).plot(jplot).fighandel_yCI) && ...
            ar.model(m).plot(jplot).fighandel_yCI ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_yCI==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_yCI, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_yCI;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.1+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_yCI = h;
        fastPlotTmp = false;
    end
else
    if(isfield(ar.model(m).plot(jplot), 'fighandel_y') && ~isempty(ar.model(m).plot(jplot).fighandel_y) && ...
            ar.model(m).plot(jplot).fighandel_y ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_y==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_y, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_y;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.1+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_y = h;
        fastPlotTmp = false;
    end
end
if(~fastPlot)
    clf
end



function savePath = mySaveFigure(h, name)
global ar
if(ar.config.ploterrors == -1)
    savePath = [arSave '/FiguresCI/Y'];
else
    savePath = [arSave '/Figures/Y'];
end

if(~exist(savePath, 'dir'))
	mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc', savePath);
% print('-dpng', savePath);
system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
% plot2svg([savePath '.svg'], h);



function str = mypath(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');



function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);



function [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols)
global ar
ny = size(ar.model(jm).data(jd).y, 2);
[nrows, ncols] = arNtoColsAndRows(ny, rowstocols);

