
% Plot models X
%
% arPlotX(saveToFile, fastPlot)
%
% saveToFile    [false]
% fastPlot      [false]

function arPlotX(saveToFile, fastPlot)

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

% constants
labelfontsize = 12;
labelfonttype = 'Arial';
rowstocols = 0.5;
overplot = 0.1;

figcount = 1;
for jm = 1:length(ar.model)
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotXs(jplot)==1)
            if(ar.config.ploterrors == -1)
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['CI-X: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            else
                [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['X: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
            end
            
            % plotting
            ccount = 1;
            if(~ar.model(jm).plot(jplot).doseresponse)
                cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                
                for jd = ar.model(jm).plot(jplot).dLink
                    [t, u, x, ulb, uub, xlb, xub, jc, dxdt] = getData(jm, jd);
                    
                    % rows and cols
                    [ncols, nrows, nu, nx, iu, ix] = myColsAndRows(jm, rowstocols);
                    
                    Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount);
                    countu = 0;
                    for ju = iu
                        countu = countu + 1;
                        if(~fastPlotTmp)
                            g = subplot(nrows,ncols,countu);
                            ar.model(jm).plot(jplot).gu(ju) = g;
                            mySubplotStyle(g, labelfontsize, labelfonttype);
                            ltmp = plot(g, t, u(:,ju), Clines{:});
                            cclegendstyles(ccount) = ltmp;
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.u(ju,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.u(ju,jc) = ltmp;
                            end
                            hold(g, 'on');
                            if(ar.config.ploterrors == -1)
                                tmpx = [t(:); flipud(t(:))];
                                tmpy = [uub(:,ju); flipud(ulb(:,ju))];
                                ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                            end
                        else
                            if(jd~=0)
                                set(ar.model(jm).data(jd).plot.u(ju,jc), 'YData', u(:,ju));
                            else
                                set(ar.model(jm).condition(jc).plot.u(ju,jc), 'YData', u(:,ju));
                            end
                        end
                    end
                    countx = 0;
                    for jx = ix
                        countx = countx + 1;
                        if(~fastPlotTmp)
                            g = subplot(nrows,ncols,countx+nu);
                            ar.model(jm).plot(jplot).gx(jx) = g;
                            mySubplotStyle(g, labelfontsize, labelfonttype);
                            
                            % plot ssa
                            if(isfield(ar.model(jm).condition(jc), 'xFineSSA'))
                                for jssa = 1:size(ar.model(jm).condition(jc).xFineSSA_lb, 3)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ar.model(jm).condition(jc).xFineSSA_ub(:,jx,jssa); ...
                                        flipud(ar.model(jm).condition(jc).xFineSSA_lb(:,jx,jssa))];
                                    patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.2+0.8, 'FaceColor', Clines{2}*0.2+0.8)
                                    %                                 patch(tmpx, tmpy, tmpx*0-2*eps, 'EdgeColor', Clines{2}*0.4+0.6, 'FaceColor', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                for jssa = 1:size(ar.model(jm).condition(jc).xFineSSA_lb, 3)
                                    plot(t, ar.model(jm).condition(jc).xFineSSA(:,jx,jssa), 'Color', Clines{2}*0.4+0.6)
                                    hold(g, 'on');
                                end
                                if(size(ar.model(jm).condition(jc).xFineSSA,3)>1)
                                    plot(t, mean(ar.model(jm).condition(jc).xFineSSA(:,jx,:),3), '--', 'Color', Clines{2})
                                end
                                
                                % density plot
                                %                             xtmp = linspace(min(min(ar.model(jm).condition(jc).xFineSSA_lb(:,jx,:)))*0.8, ...
                                %                                 max(max(ar.model(jm).condition(jc).xFineSSA_ub(:,jx,:)))*1.2, 100);
                                %                             xSSAdens = zeros(length(t), length(xtmp));
                                %                             for jssa=1:length(t)
                                %                                 xSSAdens(jssa,:) = ksdensity(squeeze((ar.model(jm).condition(jc).xFineSSA_ub(jssa,jx,:) + ...
                                %                                 ar.model(jm).condition(jc).xFineSSA_lb(jssa,jx,:))/2), xtmp);
                                %                                 xSSAdens(jssa,:) = xSSAdens(jssa,:) / sum(xSSAdens(jssa,:));
                                %                             end
                                %                             [A,B] = meshgrid(xtmp,t);
                                %                             colormap jet
                                %                             contourf(B, A, xSSAdens);
                                %                             hold(g, 'on');
                            end
                            
                            ltmp = plot(g, t, x(:,jx), Clines{:});
                            cclegendstyles(ccount) = ltmp;
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.x(jx,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.x(jx,jc) = ltmp;
                            end
                            hold(g, 'on');
                            
                            % steady state
                            xss = x(1,jx) + dxdt(jx)*(t-min(t));
                            xss(xss<0) = nan;
                            xss(xss>2*max(x(:,jx))) = nan;
                            ltmp = plot(g, t, xss, '--', Clines{:});
                            if(jd~=0)
                                ar.model(jm).data(jd).plot.xss(jx,jc) = ltmp;
                            else
                                ar.model(jm).condition(jc).plot.xss(jx,jc) = ltmp;
                            end
                            
                            if(ar.config.ploterrors == -1)
                                tmpx = [t(:); flipud(t(:))];
                                tmpy = [xub(:,jx); flipud(xlb(:,jx))];
                                ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                            end
                        else
                            % steady state
                            xss = x(1,jx) + dxdt(jx)*(t-min(t));
                            xss(xss<0) = nan;
                            xss(xss>2*max(x(:,jx))) = nan;
                            
                            if(jd~=0)
                                set(ar.model(jm).data(jd).plot.x(jx,jc), 'YData', x(:,jx));
                                set(ar.model(jm).data(jd).plot.xss(jx,jc), 'YData', xss);
                            else
                                set(ar.model(jm).condition(jc).plot.x(jx,jc), 'YData', x(:,jx));
                                set(ar.model(jm).condition(jc).plot.x(jx,jc), 'YData', xss);
                            end
                        end
                    end
                    ccount = ccount + 1;
                end
            else
                times = [];
                for jd = ar.model(jm).plot(jplot).dLink
					times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                    [ncols, nrows, nu, nx, iu, ix] = myColsAndRows(jm, rowstocols);
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
                        
                        countu = 0;
                        for ju = iu
                            countu = countu + 1;
                            [t, u, lb, ub, zero_break] = getDataDoseResponseU(jm, ju, ds, times(jt));
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,countu);
                                ar.model(jm).plot(jplot).gu(ju) = g;
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, u, Clines{:});
                                cclegendstyles(ccount) = ltmp;
                                ar.model(jm).data(jd).plot.u(ju,jt,jc) = ltmp;
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ub(:); flipud(lb(:))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                set(ar.model(jm).data(jd).plot.u(ju,jt,jc), 'YData', u);
                            end
                        end
                        countx = 0;
                        for jx = ix
                            countx = countx + 1;
                            [t, x, lb, ub, zero_break] = getDataDoseResponseX(jm, jx, ds, times(jt));
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,countx+nu);
                                ar.model(jm).plot(jplot).gx(jx) = g;
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, x, Clines{:});
                                cclegendstyles(ccount) = ltmp;
                                ar.model(jm).data(jd).plot.x(jx,jt,jc) = ltmp;
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [ub(:); flipud(lb(:))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                                if(~isempty(zero_break))
                                    plot([zero_break zero_break], ylim, 'k--');
                                end
                            else
                                set(ar.model(jm).data(jd).plot.x(jx,jt,jc), 'YData', x);
                            end
                        end
                        ccount = ccount + 1;
                    end
                end
            end
            
            % axis & titles
            for jd = ar.model(jm).plot(jplot).dLink
                countu = 0;
                for ju = iu
                    countu = countu + 1;
                    g = ar.model(jm).plot(jplot).gu(ju);
                    if(~fastPlotTmp)
                        hold(g, 'off');
                        
                        title(g, myNameTrafo(ar.model(jm).u{ju}));
                        if(ju == 1 && (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
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
                        if(nx==0 && countu == (nrows-1)*ncols + 1)
                            if(~ar.model(jm).plot(jplot).doseresponse)
                                xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                            else
                                xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                            end
                        end
                        ylabel(g, sprintf('%s [%s]', ar.model(jm).uUnits{ju,3}, ar.model(jm).uUnits{ju,2}));
                    end
                    arSpacedAxisLimits(g, overplot);
                end
                countx = 0;
                for jx = ix
                    countx = countx + 1;
                    g = ar.model(jm).plot(jplot).gx(jx);
                    if(~fastPlotTmp)
                        hold(g, 'off');
                        
                        if(nu == 0 && jx == 1 && (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
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
                        
                        if(isfield(ar.model(jm), 'xNames') && ~isempty(ar.model(jm).xNames{jx}) && ...
                                ~strcmp(ar.model(jm).xNames{jx},ar.model(jm).x{jx}))
                            title(g, [myNameTrafo(ar.model(jm).xNames{jx}) ' (' myNameTrafo(ar.model(jm).x{jx}) ')']);
                        else
                            title(g, myNameTrafo(ar.model(jm).x{jx}));                            
                        end
                        if(countx+nu == (nrows-1)*ncols + 1)
                            if(~ar.model(jm).plot(jplot).doseresponse)
                                xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                            else
                                xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                            end
                        end
                        ylabel(g, sprintf('%s [%s]', ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}));
                    end
                    arSpacedAxisLimits(g, overplot);
                end
            end
            
            if(saveToFile)
                if(ar.config.ploterrors == -1)
                    ar.model(jm).plot(jplot).savePath_FigXCI = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                else
                    ar.model(jm).plot(jplot).savePath_FigX = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                end
            end
            
            figcount = figcount + 1;
        else
            try %#ok<TRYNC>
                close(ar.model(jm).plot(jplot).fighandel_x)
            end
            ar.model(jm).plot(jplot).fighandel_x = [];
        end
    end
end

function [t, u, x, ulb, uub, xlb, xub, jc, dxdt] = getData(jm, jd)
global ar

if(jd~=0)
    jc = ar.model(jm).data(jd).cLink;
else
    jc = 1;
end
t = ar.model(jm).condition(jc).tFine;
u = ar.model(jm).condition(jc).uFineSimu;
x = ar.model(jm).condition(jc).xFineSimu;
if(isfield(ar.model(jm).condition(jc), 'xExpUB'))
    ulb = ar.model(jm).condition(jc).uFineLB;
    uub = ar.model(jm).condition(jc).uFineUB;
    xlb = ar.model(jm).condition(jc).xFineLB;
    xub = ar.model(jm).condition(jc).xFineUB;
else
    ulb = [];
    uub = [];
    xlb = [];
    xub = [];
end
dxdt = ar.model(jm).condition(jc).dxdt;
dxdt(ar.model(jm).condition(jc).qSteadyState==0) = nan;



function [t, u, lb, ub, zero_break] = getDataDoseResponseU(jm, ju, ds, ttime)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
    for jt = find(qt')
        t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(1).value)); %#ok<AGROW>
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = ds
                if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(1).value))))
                    doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(1).value)); %#ok<AGROW>
                end
            end
            doses = sort(doses);
            if(length(doses)>1)
                t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            else
                t(ccount,1) = doses(1) - 0.1; %#ok<AGROW>
                zero_break = (t(ccount,1)+doses(1))/2;
            end
        end
        u(ccount,1) = ar.model(jm).condition(jc).uExpSimu(jt,ju); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'uExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).uExpLB(jt,ju); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).uExpUB(jt,ju); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end


function [t, x, lb, ub, zero_break] = getDataDoseResponseX(jm, jx, ds, ttime)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
    for jt = find(qt')
        t(ccount,1) = log10(str2double(ar.model(jm).data(jd).condition(1).value)); %#ok<AGROW>
        if(isinf(t(ccount,1)))
            doses = [];
            for jd2 = ds
                if(~isinf(log10(str2double(ar.model(jm).data(jd2).condition(1).value))))
                    doses(end+1) = log10(str2double(ar.model(jm).data(jd2).condition(1).value)); %#ok<AGROW>
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
        x(ccount,1) = ar.model(jm).condition(jc).xExpSimu(jt,jx); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'xExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).xExpLB(jt,jx); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).xExpUB(jt,jx); %#ok<AGROW>
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
    if(isfield(ar.model(m).plot(jplot), 'fighandel_xCI') && ~isempty(ar.model(m).plot(jplot).fighandel_xCI) && ...
            ar.model(m).plot(jplot).fighandel_xCI ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_xCI==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_xCI, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_xCI;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.4+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_xCI = h;
        fastPlotTmp = false;
    end    
else
    if(isfield(ar.model(m).plot(jplot), 'fighandel_x') && ~isempty(ar.model(m).plot(jplot).fighandel_x) && ...
            ar.model(m).plot(jplot).fighandel_x ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_x==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_x, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_x;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.4+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_x = h;
        fastPlotTmp = false;
    end
end
if(~fastPlot)
    clf
end



function savePath = mySaveFigure(h, name)
global ar
if(ar.config.ploterrors == -1)
    savePath = [arSave '/FiguresCI/X'];
else
    savePath = [arSave '/Figures/X'];
end

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc2', savePath);
eval(['!ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);



function str = mypath(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');



function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);



function [ncols, nrows, nu, nx, iu, ix] = myColsAndRows(jm, rowstocols)
global ar
if(~isfield(ar.model(jm), 'qPlotU'))
    nu = size(ar.model(jm).u, 2);
    iu = 1:nu;
else
    nu = sum(ar.model(jm).qPlotU);
    iu = find(ar.model(jm).qPlotU);
end
if(~isfield(ar.model(jm), 'qPlotX'))
    nx = size(ar.model(jm).x, 2);
    ix = 1:nx;
else
    nx = sum(ar.model(jm).qPlotX);
    ix = find(ar.model(jm).qPlotX);
end
[nrows, ncols] = NtoColsAndRows(nu+nx, rowstocols);



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);


