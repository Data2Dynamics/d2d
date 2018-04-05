% Plot models Y
%
% arPlotY(saveToFile, fastPlot, doLegends, flags);
%
% saveToFile    [false]
% fastPlot      [false]
% doLegends     [true]
% flags         [none]
% 
%   After clicking the subplot of interest, the following command provides
%   annotation of the displayed plot:
%   get(gca,'UserData') 
% 
%   ar.model(jm).data(jd).highlight(jt,jy) can be used to highlight data points
%
% Flags allow you to customize the plot. HideLL, hides the -log(L) value,
% while nameOnly makes sure either the name or observable variable is used 
% as title, but not both.
% 
% Note: fastplot = 2 suppresses closing the figure and plots a new figure 
% for each plot.
% 
% 
% Doku: 
% https://github.com/Data2Dynamics/d2d/wiki/Plotting-options-and-the-meaning-of-ar.config.ploterrors

function arPlotY(saveToFile, fastPlot, doLegends, varargin)

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
switches = {'hidell', 'nameonly'};
[opts] = argSwitch( switches, varargin );
if(isfield(ar.config,'nfine_dr_plot'))
    nfine_dr_plot = ar.config.nfine_dr_plot;
    nfine_dr_method = ar.config.nfine_dr_method;
else
    nfine_dr_plot = 1;
    nfine_dr_method = 'spline';
end
% if(~isfield(ar.config,'useFitErrorMatrix'))
%     ar.config.useFitErrorMatrix = false;
% end

figcount = 1;
% Set ar.config.debugExp to 1 to mark points on the simulation curve with
% the values that are actually used in the likelihood calculation.
expPlot = isfield(ar.config, 'debugExp')&&ar.config.debugExp;
for jm = 1:length(ar.model)
    ar.model(jm).chi2 = 0;
    ar.model(jm).ndata = 0;
    
    for jplot = 1:length(ar.model(jm).plot)
        isBarGraph = 0;
        if(ar.model(jm).qPlotYs(jplot)==1 && ar.model(jm).plot(jplot).ny>0)
%             if( (ar.config.useFitErrorMatrix == 0 && ar.config.ploterrors == -1) || ...
%                     (ar.config.useFitErrorMatrix == 1 && ar.config.ploterrors_matrix(jm,ar.model(jm).plot(jplot).dLink(1))==-1) )
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['Y: ploterrors' num2str(ar.config.ploterrors) ar.model(jm).plot(jplot).name], figcount, fastPlot);
%             else
%                 [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['Y: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
%             end
            
            if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
                logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
            else
                logplotting_xaxis = true;
            end
            
            % plotting
            ccount = 1;
            chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
            ndata = zeros(1,ar.model(jm).plot(jplot).ny);   
            if(~ar.model(jm).plot(jplot).doseresponse)
              %% Time courses
                cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                
                for jd = ar.model(jm).plot(jplot).dLink
                    % rows and cols
                    ny = size(ar.model(jm).data(jd).y, 2);
                    [nrows, ncols] = arNtoColsAndRows(ny);
                    if(nrows*ncols == ny)
                        [nrows, ncols] = arNtoColsAndRows(ny+1);
                    end
                    
                    Clines = arLineMarkersAndColors(ccount, ...
                        length(ar.model(jm).plot(jplot).dLink), ...
                        [], 'none', '-');
                    ClinesExp = arLineMarkersAndColors(ccount, ...
                        length(ar.model(jm).plot(jplot).dLink), ...
                        [], 'none', 'none');
                    
                    for jy = 1:ny
                        whichYplot = arWhichYplot(jm,jd,[],jy);
                        [t, y, ystd, tExp, yExp, yExpStd, lb, ub, yExpHl, yExpSimu] = getData(jm, jd, jy);
                        
                        % Display bar for single time point measurements
                        if ( isfield( ar.config, 'barhack' ) && ( ar.config.barhack == 1 ) )
                            if ( numel( unique( t(~isnan(t)) ) ) == 1 )
                                t = t(~isnan(t));
                                y = nanmean(y);
                                ystd = nanmean(ystd);
                                t = [t-0.01, t-0.01, t+0.01, t+0.01];
                                y = [0; y; y; 0];
                                ystd = [0; ystd; ystd; 0];
                                isBarGraph = 1;
                                set(gca, 'YGrid', 'on');
                                set(gca, 'XTick', []);
                            end
                        end
                                           
                        if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                            trafo = @(x)10.^x;
                        else
                            trafo = @(x)x;
                        end
                        if(~fastPlotTmp)
                            g = subplot(nrows,ncols,jy);
                            ar.model(jm).plot(jplot).gy(jy) = g;
                            plotResFuncSpecificElements( g, t, jm, jd, jy, trafo, Clines{2} );
                            
                            if(~isfield(ar.model(jm).data(jd),'qFit') || ar.model(jm).data(jd).qFit(jy))
                                ClinesExp{6} = '*';
                            else
                                ClinesExp{6} = 'o';
                            end
                                
                            % plot ssa
                            if isfield(ar.model(jm).data(jd), 'yFineSSA') % && ar.config.useFitErrorMatrix == 0 && ar.config.ploterrors==1) || ...
%                                         (isfield(ar.model(jm).data(jd), 'yFineSSA') && ar.config.useFitErrorMatrix==1 && ar.config.ploterrors_matrix(jm,jd)==1) )
                                if(isfield(ar.model(jm).data(jd), 'yFineSSA_lb'))
                                    for jssa = 1:size(ar.model(jm).data(jd).yFineSSA_lb, 3)
                                        tmpx = [t(:); flipud(t(:))];
                                        tmpy = trafo([ar.model(jm).data(jd).yFineSSA_ub(:,jy,jssa); ...
                                            flipud(ar.model(jm).data(jd).yFineSSA_lb(:,jy,jssa))]);
                                        patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
                                        hold(g, 'on');

                                    end
                                end
                                for jssa = 1:size(ar.model(jm).data(jd).yFineSSA, 3)
                                    plot(t, trafo(ar.model(jm).data(jd).yFineSSA(:,jy,jssa)), 'Color', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                if(size(ar.model(jm).data(jd).yFineSSA,3)>1)
                                    plot(t, trafo(mean(ar.model(jm).data(jd).yFineSSA(:,jy,:),3)), '--', 'Color', Clines{2})
                                end
                            end

                            tmpx = t;
                            tmpy = trafo(y);
                            qfinite = ~isinf(tmpy);
                            if(sum(qfinite)>0)
                                ar.model(jm).data(jd).plot.y(jy) = plot(g, tmpx(qfinite), tmpy(qfinite), Clines{:});
                                cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                            end
                            hold(g, 'on');
                            if any(whichYplot==3:5)
                                tmpx = [t(:); flipud(t(:))];
                                if any(whichYplot==[3,5])
                                    tmpy = trafo([y + ystd; flipud(y - ystd)]);
                                elseif whichYplot==4 && ~isempty(ub)
                                    tmpy = trafo([ub; flipud(lb)]);
                                else
                                    tmpy = NaN(size(tmpx));
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
                                if any(whichYplot==[1,3,4])
                                    plot(g, tExp, trafo(yExp), ClinesExp{:});
                                    if(sum(~isnan(yExpHl))>0)
                                        hold(g,'on');
                                        plot(g, tExp, trafo(yExpHl), ClinesExp{:},'LineWidth',2,'MarkerSize',10);
                                    end
                                else
                                    errorbar(g, tExp, trafo(yExp), trafo(yExp) - trafo(yExp - yExpStd), trafo(yExp + yExpStd) - trafo(yExp), ClinesExp{:});
                                    if(sum(~isnan(yExpHl))>0)
                                        hold(g,'on');
                                        errorbar(g, tExp, trafo(yExpHl), trafo(yExp) - trafo(yExp - yExpStd), trafo(yExp + yExpStd) - trafo(yExp), ClinesExp{:},'LineWidth',2,'MarkerSize',10);
                                    end
                                end
                                if ( expPlot )
                                    plot(g, tExp, trafo(yExpSimu), 'x', ClinesExp{1:2}, 'Markersize', 6);
                                end
                            end
                        else
                            tmpy = trafo(y);
                            qfinite = ~isinf(tmpy);
                            set(ar.model(jm).data(jd).plot.y(jy), 'YData', tmpy(qfinite));
                            if any(whichYplot==[3,5])
                                tmpy = trafo([y + ystd; flipud(y - ystd)]);
                                qfinite = ~isinf(tmpy);
                                if(sum(qfinite)>0)
                                    set(ar.model(jm).data(jd).plot.ystd(jy),  'YData', tmpy(qfinite));
                                    set(ar.model(jm).data(jd).plot.ystd2(jy), 'YData', tmpy(qfinite));
                                end
                            end
                        end
                        
                        % chi^2 & ndata
                        if(isfield(ar.model(jm).data(jd),'qFit') && ar.model(jm).data(jd).qFit(jy)==1 && ~isempty(ar.model(jm).data(jd).chi2))
                            chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                            ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                            if ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0) % ( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1) || ...
%                                     (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                                chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2err(jy);
                            end
                        end
                    end
                    ccount = ccount + 1;
                end
            else
              %% Dose responses
                times = [];
                for jd = ar.model(jm).plot(jplot).dLink
                    times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                    ny = size(ar.model(jm).data(jd).y, 2);
                    [nrows, ncols] = arNtoColsAndRows(ny);
                    if(nrows*ncols == ny)
                        [nrows, ncols] = arNtoColsAndRows(ny+1);
                    end

                    for jy = 1:ny
                        % chi^2 & ndata
                        if(ar.model(jm).data(jd).qFit(jy)==1)
                            chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                            ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                            if ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0)
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
                        Clines = arLineMarkersAndColors(ccount, ...
                            length(times)*length(jcs), ...
                            [], 'none', '-');
                        ClinesExp = arLineMarkersAndColors(ccount, ...
                            length(times)*length(jcs), ...
                            [], 'none', 'none');
                        
                        for jy = 1:ny
                            whichYplot = arWhichYplot(jm,jd,[],jy);
                            
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, data_qFit, yExpHl] = ...
                                getDataDoseResponse(jm, jy, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
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
                                if(~isempty(lb))
                                    lb = lb(qit);
                                    ub = ub(qit);
                                    lb = interp1(t,lb,tf,nfine_dr_method);
                                    ub = interp1(t,ub,tf,nfine_dr_method);
                                end
                                t = tf;
                            end
                           
                            if any(whichYplot==[3]) % ( (ar.config.useFitErrorMatrix == 0 && ar.config.ploterrors==0) || ...
%                                     (ar.config.useFitErrorMatrix==1 && ar.config.ploterrors_matrix(jm,jd)==0) )
                                lb = y(:) - ystd(:);
                                ub = y(:) + ystd(:);
                            end

                            if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                trafo = @(x)10.^x;
                            else
                                trafo = @(x)x;
                            end
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,jy);
                                ar.model(jm).plot(jplot).gy(jy) = g;
                                plotResFuncSpecificElements( g, t, jm, jd, jy, trafo, Clines{2} );
                                
                                if(data_qFit)
                                    ClinesExp{6} = '*';
                                else
                                    ClinesExp{6} = 'o';
                                end
                                
                                qfinite = ~isinf(t) & ~isinf(y);
                                ar.model(jm).data(jd).plot.y(jy,jt,jc) = plot(g, t(qfinite), trafo(y(qfinite)), Clines{:});
                                cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy,jt,jc);
                                hold(g, 'on');

                                if any(whichYplot==[3,5]) % 5?
                                   tmpx = [t(:); flipud(t(:))];
                                   tmpy = trafo([ub; flipud(lb)]);
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
                                    if any(whichYplot==[1,3,4])
                                        plot(g, tExp, trafo(yExp), ClinesExp{:});
                                        if(sum(~isnan(yExpHl))>0)
                                            hold(g,'on');
                                            plot(g, tExp, trafo(yExpHl), ClinesExp{:}, 'LineWidth',2, 'MarkerSize', 10);                                            
                                        end
                                    else
                                        errorbar(g, tExp, trafo(yExp), trafo(yExp) - trafo(yExp - yExpStd), trafo(yExp + yExpStd) - trafo(yExp), ClinesExp{:});
                                        if(sum(~isnan(yExpHl))>0)
                                            hold(g,'on')
                                            errorbar(g, tExp, trafo(yExpHl), trafo(yExp) - trafo(yExp - yExpStd), trafo(yExp + yExpStd) - trafo(yExp), ClinesExp{:}, 'LineWidth', 2, 'MarkerSize', 10);
                                        end
                                    end
                                end
                                
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                ytmp = trafo(y);
								qfinite = ~isinf(t) & ~isinf(ytmp);
                                set(ar.model(jm).data(jd).plot.y(jy,jt,jc), 'YData', ytmp(qfinite));
                                if any(whichYplot==4) && ~isempty(ub)
                                    tmpx = [t(:); flipud(t(:))];
									tmpy = trafo([ub; flipud(lb)]);
                                    qfinite = ~isinf(tmpy) & ~isinf(tmpx);
                                    if(sum(qfinite)>0)
                                        set(ar.model(jm).data(jd).plot.ystd(jy,jt,jc),  'YData', tmpy(qfinite));
                                        set(ar.model(jm).data(jd).plot.ystd2(jy,jt,jc), 'YData', tmpy(qfinite));
                                    end
                                end
                            end
                        end
                        ccount = ccount + 1;
                    end
                end
            end
            
            if ( isBarGraph )
                set(gca, 'YGrid', 'on');
                set(gca, 'XTick', []);                
            end
            
            % axis & titles
            if(~fastPlotTmp && exist('suptitle','file')==2 && isfield(ar.config, 'useSuptitle') && ar.config.useSuptitle) % suptitle function is available
                suptitle(arNameTrafo([ar.model(jm).name,': ',ar.model(jm).plot(jplot).name]))
            end
            
            for jc = 1:length(ar.model(jm).data(jd).condition)
                if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
                    jcondi = jc;
                end
            end
            for jy = 1:ny
                g = ar.model(jm).plot(jplot).gy(jy);
                if(~fastPlotTmp)
                    hold(g, 'off');
                    arSubplotStyle(g);
                    
                    if ( ~isBarGraph )
                        qxlabel = jy == (nrows-1)*ncols + 1;
                        if(ny <= (nrows-1)*ncols)
                            qxlabel = jy == (nrows-2)*ncols + 1;
                        end
                        if(qxlabel)
                            if(~ar.model(jm).plot(jplot).doseresponse)
                                xlabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).tUnits{3}, ar.model(jm).data(jd).tUnits{2}));
                            else
                                if(isfield(ar.model(jm).plot(jplot), 'response_parameter') && ...
                                        ~isempty(ar.model(jm).plot(jplot).response_parameter))
                                    resppar = ar.model(jm).plot(jplot).response_parameter;
                                else
                                    resppar = arNameTrafo(ar.model(jm).data(jd).condition(jcondi).parameter);
                                end
                                if(logplotting_xaxis)
                                    xlabel(g, sprintf('log_{10}(%s)', resppar));
                                else
                                    xlabel(g, sprintf('%s', resppar));
                                end
                            end
                        end
                    end
                    if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
                        ylabel(g, sprintf('log_{10}(%s) [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    else
                        ylabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    end
                    
                    if(doLegends && jy == ny && (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                        hl = [];
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            if(length(ar.model(jm).plot(jplot).dLink)>1)
                                hl = legend(g, cclegendstyles, arNameTrafo(ar.model(jm).plot(jplot).condition));
                            end
                        else
                            legendtmp = {};
                            ccount = 1;
                            for jt=1:length(times)
                                if(~isempty(conditions))
                                    for jc = 1:length(conditions)
                                        if(~isempty(conditions{jc}))
                                            legendtmp{ccount} = sprintf('t=%g%s : %s', times(jt), ar.model(jm).tUnits{2}, conditions{jc}); %#ok<AGROW>
                                        else
                                            legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                        end
                                        ccount = ccount + 1;
                                    end
                                else
                                    legendtmp{ccount} = sprintf('t=%g%s', times(jt), ar.model(jm).tUnits{2}); %#ok<AGROW>
                                    ccount = ccount + 1;
                                end
                            end
                            hl = legend(g, cclegendstyles, arNameTrafo(legendtmp));
                        end
                        if(~isempty(hl))
                            gref = subplot(nrows,ncols,nrows*ncols);
                            axis(gref,'off');
                            grefloc = get(gref, 'Position');
                            legloc = get(hl, 'Position');
                            legloc(1:2) = grefloc(1:2) + [.95*(grefloc(3)-legloc(3)) grefloc(4)-legloc(4)];
                            % legloc(1:2) = [1-legloc(3) 0];
                            if(sum(isinf(legloc))==0)
                                set(hl, 'Position', legloc);
                            end
                        end
                    end
                end
                titstr = {};
                if(isfield(ar.model(jm).data(jd), 'yNames') && ~isempty(ar.model(jm).data(jd).yNames{jy}) && ...
                        ~strcmp(ar.model(jm).data(jd).yNames{jy}, ar.model(jm).data(jd).y{jy}))
                    if ( opts.nameonly)
                        titstr{1} = arNameTrafo(ar.model(jm).data(jd).yNames{jy});
                    else
                        titstr{1} = [arNameTrafo(ar.model(jm).data(jd).yNames{jy}) ' (' arNameTrafo(ar.model(jm).data(jd).y{jy}) ')'];
                    end
                else
                    titstr{1} = arNameTrafo(ar.model(jm).data(jd).y{jy});
                end
                if(~opts.hidell)
                    if(isfield(ar.model(jm).data(jd), 'yExp'))
                        if(ndata(jy)>0)
                            if any(whichYplot==[2,5]) % ( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1) || ...
%                                     (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                               titstr{2} = sprintf('-2 log(L)_{%i} = %g', ndata(jy), 2*ndata(jy)*log(sqrt(2*pi)) + chi2(jy));
                            else
                               titstr{2} = sprintf('\\chi^2_{%i} = %g', ndata(jy), chi2(jy));
                            end
                       end
                    end
                end
                title(g, titstr);
                set(g,'UserData',...
                    struct('jm',jm,'jplot',jplot,'jy',jy, ...
                    'dLink',ar.model(jm).plot(jplot).dLink, ...
                    'yName',ar.model(jm).data(jd).yNames{jy}, ...
                    'model_name',ar.model(jm).name, ...
                    'plot_name',ar.model(jm).plot(jplot).name, ...
                    'data_name',ar.model(jm).data(jd).name, ...
                    'fy',ar.model(jm).data(jd).fy{jy}, ...
                    'fystd',ar.model(jm).data(jd).fystd{jy} ...                    
                    ))
                arSpacedAxisLimits(g);
            end
            
            ar.model(jm).plot(jplot).chi2 = sum(chi2);
            ar.model(jm).plot(jplot).ndata = sum(ndata);
            
            ar.model(jm).chi2 = ar.model(jm).chi2 + sum(chi2);
            ar.model(jm).ndata = ar.model(jm).ndata + sum(ndata);
            
            figcount = figcount + 1;
            
%             if ( ~isfield( ar.config, 'useFitErrorMatrix' ) )
%                 ar.config.useFitErrorMatrix = 0;
%             end
            
            if(saveToFile)                
%                 if any(whichYplot==[2,4,5])% ( (ar.config.useFitErrorMatrix==0 && ar.config.ploterrors == -1) || ...
% %                         (ar.config.useFitErrorMatrix==1 && ar.config.ploterrors_matrix(jm,jd)==-1) )                        
                    [ ar.model(jm).plot(jplot).savePath_FigY , ...
                      ar.model(jm).plot(jplot).nRows, ...
                      ar.model(jm).plot(jplot).nCols ] = ...
                        arSaveFigure(h, ar.model(jm).plot(jplot).name, ['/Figures/Yold_ploterrors' num2str(ar.config.ploterrors) 'fiterrors' num2str(ar.config.fiterrors)]);
%                 else
%                     [ ar.model(jm).plot(jplot).savePath_FigY, ...
%                       ar.model(jm).plot(jplot).nRows, ...
%                       ar.model(jm).plot(jplot).nCols ] = ...
%                          arSaveFigure(h, ar.model(jm).plot(jplot).name, '/Figures/Y');
%                 end
            end
        else
            try %#ok<TRYNC>
                if (fastplot ~= 2)
                    close(ar.model(jm).plot(jplot).fighandel_y)
                end
            end
            ar.model(jm).plot(jplot).fighandel_y = [];
        end
    end
    
    
end

function plotResFuncSpecificElements( g, t, jm, jd, jy, trafo, col )
    global ar;
    if( isfield(ar.model(jm).data(jd), 'resfunction') )
        if ( isstruct( ar.model(jm).data(jd).resfunction ) )
            if ( strcmp( ar.model(jm).data(jd).resfunction.type, 'DetectionLimit' ) )
                plot(g, [min(t), max(t)], trafo([ar.model(jm).data(jd).resfunction.LoD(jy), ar.model(jm).data(jd).resfunction.LoD(jy)]), ':', 'Color', min( col + [.2, .2, .2], [1, 1, 1]));
                hold(g, 'on' );
            end
        end
    end

function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, yExpHl, yExpSimu] = getData(jm, jd, jy)
global ar

if(isfield(ar.model(jm).data(jd),'tFine'))
    t = ar.model(jm).data(jd).tFine;
    y = ar.model(jm).data(jd).yFineSimu(:,jy);
    ystd = ar.model(jm).data(jd).ystdFineSimu(:,jy);
else
    t = nan;
    y = nan;
    ystd = nan;
end
if(isfield(ar.model(jm).data(jd), 'yExp') && ~isempty(ar.model(jm).data(jd).yExp))
    tExp = ar.model(jm).data(jd).tExp;
    yExp = ar.model(jm).data(jd).yExp(:,jy);
    yExpSimu = ar.model(jm).data(jd).yExpSimu(:,jy);
    if ar.config.fiterrors==-1
        yExpStd = ar.model(jm).data(jd).yExpStd(:,jy);
    elseif any(ar.config.fiterrors == [0,1])
        yExpStd = nan*ar.model(jm).data(jd).yExpStd(:,jy);
        
        if  ar.config.ploterrors==1 || ar.config.fiterrors==1% error model as error bars only if "ploterrors==1" 
            if(isfield(ar.model(jm).data(jd),'ystdExpSimu'))
                yExpStd = ar.model(jm).data(jd).ystdExpSimu(:,jy);
            end
        end
        
        if ar.config.fiterrors == 0
            notnan = ~isnan(ar.model(jm).data(jd).yExpStd(:,jy));
            yExpStd(notnan) = ar.model(jm).data(jd).yExpStd(notnan,jy);
        end
    end
    if(isfield(ar.model(jm).data(jd),'highlight'))
        hl = ar.model(jm).data(jd).highlight(:,jy);
    else
        hl = zeros(size(yExp));
    end
    yExpHl = yExp;
    yExpHl(hl==0) = NaN;
else
    tExp = [];
    yExp = [];
	yExpStd = [];
    yExpHl = [];
end
if(isfield(ar.model(jm).data(jd), 'yFineLB'))
	lb = ar.model(jm).data(jd).yFineLB(:,jy);
	ub = ar.model(jm).data(jd).yFineUB(:,jy);
else
	lb = [];
	ub = [];
end
if(sum(abs(size(yExp)-size(yExpHl)))>0)
    error('')
end



function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, data_qFit, yExpHl] = ...
    getDataDoseResponse(jm, jy, ds, ttime, dLink, logplotting_xaxis)
global ar


zero_break = [];
data_qFit = true;

ccount = 1;
for jd = ds
    data_qFit = data_qFit & ar.model(jm).data(jd).qFit(jy);
	qt = ar.model(jm).data(jd).tExp == ttime;
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    for jt = find(qt')
        if(logplotting_xaxis)
            t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(jcondi).value)); %#ok<AGROW>
        else
            t(ccount,1) = str2double(ar.model(jm).data(jd).condition(jcondi).value); %#ok<AGROW>
        end
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = dLink
                if(logplotting_xaxis)
                    if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value))))
                        doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(jcondi).value)); %#ok<AGROW>
                    end
                else
                    doses(end+1) = str2double(ar.model(jm).data(jd2).condition(jcondi).value); %#ok<AGROW>
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
        yExpHl(ccount,1) = NaN;
        if(isfield(ar.model(jm).data(jd),'highlight'))
            if(ar.model(jm).data(jd).highlight(jt,jy)~=0)
                yExpHl(ccount,1) = yExp(ccount,1);                
            end
        end
        
        if ar.config.fiterrors==-1
            yExpStd(ccount,1) = ar.model(jm).data(jd).yExpStd(jt,jy);
        elseif any(ar.config.fiterrors == [0,1])
            yExpStd(ccount,1) = nan*ar.model(jm).data(jd).yExpStd(jt,jy);
            
            if  ar.config.ploterrors==1 || ar.config.fiterrors==1% error model as error bars only if "ploterrors==1"
                if(isfield(ar.model(jm).data(jd),'ystdExpSimu'))
                    yExpStd(ccount,1) = ar.model(jm).data(jd).ystdExpSimu(jt,jy);
                end
            end
            
            if ar.config.fiterrors == 0
                if ~isnan(ar.model(jm).data(jd).yExpStd(jt,jy));
                    yExpStd(ccount,1) = ar.model(jm).data(jd).yExpStd(jt,jy);
                end
            end
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

[tExp,itexp] = sort(tExp);
yExp = yExp(itexp);
yExpHl = yExpHl(itexp);
yExpStd = yExpStd(itexp);

[t,it] = sort(t);
y = y(it);
ystd = ystd(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end

%% sub-functions


function [h, fastPlotTmp] = myRaiseFigure(m, jplot, figname, figcount, fastPlot)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;
fastPlotTmp = fastPlot;

if ar.config.fiterrors==-1 % ( (ar.config.useFitErrorMatrix == 0 && ar.config.ploterrors == -1) || ...
%         (ar.config.useFitErrorMatrix == 1 && ar.config.ploterrors_matrix(m,ar.model(m).plot(jplot).dLink(1))==-1) )
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
if(fastPlot==0)
    clf
end

function [opts] = argSwitch( switches, varargin )

    for a = 1 : length(switches)
        opts.(lower(switches{a})) = 0;
    end
    
    if ~isempty( varargin{1} )
        for a = 1 : length( varargin{1} )
            if ( max( strcmp( lower(switches), lower(varargin{1}{a}) ) ) == 0 )
                error( 'Invalid switch argument was provided %s', varargin{1}{a} );
            end
            opts.(lower(varargin{1}{a})) = 1;
        end
    end
