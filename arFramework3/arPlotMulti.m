% Plot models and datasets
%
% arPlotMulti(ps, ps_weigths, saveToFile, filenameAddition)
%
% ps                    par values
% ps_weigths            par weights
% saveToFile            [false]
% filenameAddition      ['']

function arPlotMulti(ps, ps_weigths, saveToFile, filenameAddition)

matVer = ver('MATLAB');

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('saveToFile','var'))
	saveToFile = false;
end
if(~exist('filenameAddition','var'))
	filenameAddition = '';
end
if(~exist('ps_weigths','var') || isempty(ps_weigths))
	ps_weigths = ones(1,size(ps,1));
end

% constants
labelfontsize = 12;
labelfonttype = 'TimesNewRoman';
rowstocols = 0.5; %0.7; 0.45;
overplot = 0.1;
linesize = 0.5;

logplotting_xaxis = true;

figcount = 1;
figcountx = 1;
figcountv = 1;

% for reseting
pReset = ar.p;

np = size(ps,1);
jeti = jet(np);

for jp=1:np
    if(jp==2)
        hbar = waitbar(0, 'Please wait...');
    elseif(jp>2)
        hbar = waitbar(jp / np, hbar, 'Please wait...');
    end
    
    ar.p = ps(jp,:) + 0;
%     try
% 		arFitObs(true);
        arSimu(false, true);
        arSimu(false, false);
        weight = ps_weigths(jp);
        
        for jm = 1:length(ar.model)
            for jplot = 1:length(ar.model(jm).plot)
                if(ar.model(jm).qPlotYs(jplot) && ar.model(jm).plot(jplot).ny>0)
                    if(jp==1)
                        myRaiseFigure(jm, jplot, ['Y: ' ar.model(jm).plot(jplot).name], figcount);
                    end
                    
                    % plotting
                    ccount = 1;
                    if(~ar.model(jm).plot(jplot).doseresponse)
                        cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                        
                        for jd = ar.model(jm).plot(jplot).dLink
                            % rows and cols
                            [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols);
                            ar.model(jm).plot(jplot).ny = ny;
                            
                            Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount, weight);
                            if(length(ar.model(jm).plot(jplot).dLink)==1)
                                Clines{2} = jeti(jp,:);
                            end
                            
                            for jy = 1:ny
                                [t, y, ystd, tExp, yExp, yExpStd] = getData(jm, jd, jy); %#ok<ASGLU>
                                
                                if(jp==1)
                                    g = subplot(nrows,ncols,jy);
                                    ar.model(jm).plot(jplot).gy(jy) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gy(jy);
                                end
                                
                                if(ar.model(jm).data(jd).qFit(jy))
                                    markerstyle = '*';
                                else
                                    markerstyle = 'o';
                                end
                                
                                Clinesdata = Clines;
                                if(length(ar.model(jm).plot(jplot).dLink)==1)
                                    Clinesdata{2} = [0 0 0];
                                end
                                if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                    ar.model(jm).data(jd).plot.y(jy) = plot(g, t, 10.^y, Clines{:}, 'LineWidth', linesize);
                                    cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                    hold(g, 'on');
                                    if(isfield(ar.model(jm).data(jd), 'yExp') && jp==np)
                                        
                                        if(ar.config.ploterrors==0)
                                            plot(g, tExp, 10.^yExp, markerstyle, Clinesdata{:});
                                        else
                                            errorbar(g, tExp, 10.^yExp, 10.^yExp - 10.^(yExp - yExpStd), ...
                                                10.^(yExp + yExpStd) - 10.^yExp, markerstyle, Clinesdata{:});
                                        end
                                    end
                                else
                                    tmpx = t;
                                    tmpy = y;
                                    qfinite = ~isinf(tmpy);
                                    ar.model(jm).data(jd).plot.y(jy) = plot(g, tmpx(qfinite), tmpy(qfinite), Clines{:}, 'LineWidth', linesize);
                                    cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                    hold(g, 'on');
                                    if(isfield(ar.model(jm).data(jd), 'yExp') && jp==np)
                                        if(ar.config.ploterrors==0)
                                            plot(g, tExp, yExp, markerstyle, Clinesdata{:});
                                        else
                                            errorbar(g, tExp, yExp, yExpStd, markerstyle, Clinesdata{:});
                                        end
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
                                Clines = myLineStyle(length(times)*length(jcs), ccount, weight);
                                if(isempty(ar.model(jm).plot(jplot).condition))
                                    Clines{2} = jeti(jp,:);
                                end
                                
                                for jy = 1:ny
                                    [t, y, ~, tExp, yExp, yExpStd, ~, ~, zero_break] = ...
                                        getDataDoseResponse(jm, jy, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                                    
                                    if(jp==1)
                                        g = subplot(nrows,ncols,jy);
                                        ar.model(jm).plot(jplot).gy(jy) = g;
                                    else
                                        g = ar.model(jm).plot(jplot).gy(jy);
                                    end
                                    
                                    if(ar.model(jm).data(jd).qFit(jy))
                                        markerstyle = '*';
                                    else
                                        markerstyle = 'o';
                                    end
                                    
                                    Clinesdata = Clines;
                                    if(isempty(ar.model(jm).plot(jplot).condition))
                                        Clinesdata{2} = [0 0 0];
                                    end
                                    if(ar.model(jm).data(jd).logfitting(jy) && ~ar.model(jm).data(jd).logplotting(jy))
                                        ar.model(jm).data(jd).plot.y(jy) = plot(g, t, 10.^y, Clines{:}, 'LineWidth', linesize);
                                        cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                        hold(g, 'on');
                                        if(isfield(ar.model(jm).data(jd), 'yExp') && jp==np)
                                            if(ar.config.ploterrors==0)
                                                plot(g, tExp, 10.^yExp, markerstyle, Clinesdata{:});
                                            else
                                                errorbar(g, ar.model(jm).data(jd).tExp, 10.^yExp, ...
                                                    10.^yExp - 10.^(yExp - yExpStd), 10.^(yExp + yExpStd) - 10.^yExp, markerstyle, Clinesdata{:});
                                            end
                                        end
                                    else
                                        tmpx = t;
                                        tmpy = y;
                                        qfinite = ~isinf(tmpy);
                                        ar.model(jm).data(jd).plot.y(jy) = plot(g, tmpx(qfinite), tmpy(qfinite), Clines{:}, 'LineWidth', linesize);
                                        cclegendstyles(ccount) = ar.model(jm).data(jd).plot.y(jy);
                                        hold(g, 'on');
                                        if(isfield(ar.model(jm).data(jd), 'yExp') && jp==np)
                                            if(ar.config.ploterrors==0)
                                                plot(g, tExp, yExp, markerstyle, Clinesdata{:});
                                            else
                                                errorbar(g, tExp, yExp, yExpStd, markerstyle, Clinesdata{:});
                                            end
                                        end
                                    end
                                    if(jp==np && ~isempty(zero_break))
                                        plot(g, [zero_break zero_break], ylim(g), 'k--');
                                    end
                                end
                                ccount = ccount + 1;
                            end
                        end
                    end
                    
                    if(jp==1)
                        figcount = figcount + 1;
                    end
                end
                
                if(ar.model(jm).qPlotXs(jplot))
                    if(jp==1)
                        myRaiseFigureX(jm, jplot, ['X: ' ar.model(jm).plot(jplot).name], figcountx);
                    end
                    
                    % plotting
                    ccount = 1;
                    if(~ar.model(jm).plot(jplot).doseresponse)
                        cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                        
                        for jd = ar.model(jm).plot(jplot).dLink
                            [t, u, x, jc] = getDataX(jm, jd);
                            
                            % rows and cols
                            [ncols, nrows, nu, ~, iu, ix] = myColsAndRowsX(jm, rowstocols);
                            
                            Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount, weight);
                            if(length(ar.model(jm).plot(jplot).dLink)==1)
                                Clines{2} = jeti(jp,:);
                            end
                            
                            countu = 0;
                            for ju = iu
                                countu = countu + 1;
                                if(jp==1)
                                    g = subplot(nrows,ncols,countu);
                                    ar.model(jm).plot(jplot).gu(ju) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gu(ju);
                                end
                                
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, u(:,ju), Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                if(jd~=0)
                                    ar.model(jm).data(jd).plot.u(ju,jc) = ltmp;
                                else
                                    ar.model(jm).condition(jc).plot.u(ju,jc) = ltmp;
                                end
                                hold(g, 'on');
                            end
                            countx = 0;
                            for jx = ix
                                countx = countx + 1;
                                if(jp==1)
                                    g = subplot(nrows,ncols,countx+nu);
                                    ar.model(jm).plot(jplot).gx(jx) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gx(jx);
                                end
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, x(:,jx), Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                if(jd~=0)
                                    ar.model(jm).data(jd).plot.x(jx,jc) = ltmp;
                                else
                                    ar.model(jm).condition(jc).plot.x(jx,jc) = ltmp;
                                end
                                hold(g, 'on');
                            end
                            ccount = ccount + 1;
                        end
                    else
                        times = [];
                        ds = ar.model(jm).plot(jplot).dLink;
                        for jd = ds
							times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                            [ncols, nrows, nu, ~, iu, ix] = myColsAndRowsX(jm, rowstocols);
                        end
                        
                        cclegendstyles = zeros(1,length(times));
                        for jt = 1:length(times)
                            Clines = myLineStyle(length(times), jt, weight);
                            if(isempty(ar.model(jm).plot(jplot).condition))
                                Clines{2} = jeti(jp,:);
                            end
                            
                            countu = 0;
                            for ju = iu
                                countu = countu + 1;
                                [t, u, ~, ~, zero_break] = getDataDoseResponseU(jm, ju, ds, times(jt));
                                if(jp==1)
                                    g = subplot(nrows,ncols,countu);
                                    ar.model(jm).plot(jplot).gu(ju) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gu(ju);
                                end
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, u, Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                ar.model(jm).data(jd).plot.u(ju,jt) = ltmp;
                                hold(g, 'on');
                                if(jp==np && ~isempty(zero_break))
                                    plot(g, [zero_break zero_break], ylim(g), 'k--');
                                end
                            end
                            countx = 0;
                            for jx = ix
                                countx = countx + 1;
                                [t, x, ~, ~, zero_break] = getDataDoseResponseX(jm, jx, ds, times(jt));
                                if(jp==1)
                                    g = subplot(nrows,ncols,countx+nu);
                                    ar.model(jm).plot(jplot).gx(jx) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gx(jx);
                                end
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, x, Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                ar.model(jm).data(jd).plot.x(jx,jt) = ltmp;
                                hold(g, 'on');
                                if(jp==np && ~isempty(zero_break))
                                    plot(g, [zero_break zero_break], ylim(g), 'k--');
                                end
                            end
                        end
                    end
                    
                    if(jp==1)
                        figcountx = figcountx + 1;
                    end
                end
                
                if(ar.model(jm).qPlotVs(jplot))
                    if(jp==1)
                        myRaiseFigureV(jm, jplot, ['V: ' ar.model(jm).plot(jplot).name], figcountv);
                    end
                    
                    % plotting
                    ccount = 1;
                    if(~ar.model(jm).plot(jplot).doseresponse)
                        cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                        
                        for jd = ar.model(jm).plot(jplot).dLink
                            [t, v, jc] = getDataV(jm, jd);
                            
                            % rows and cols
                            [ncols, nrows, iv] = myColsAndRowsV(jm, rowstocols);
                            
                            Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount, weight);
                            if(length(ar.model(jm).plot(jplot).dLink)==1)
                                Clines{2} = jeti(jp,:);
                            end
                            
                            countv = 0;
                            for jv = iv
                                countv = countv + 1;
                                if(jp==1)
                                    g = subplot(nrows,ncols,countv);
                                    ar.model(jm).plot(jplot).gv(jv) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gv(jv);
                                end
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, v(:,jv), Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                if(jd~=0)
                                    ar.model(jm).data(jd).plot.v(jv,jc) = ltmp;
                                else
                                    ar.model(jm).condition(jc).plot.v(jv,jc) = ltmp;
                                end
                                hold(g, 'on');
                            end
                            ccount = ccount + 1;
                        end
                    else
                        times = [];
                        ds = ar.model(jm).plot(jplot).dLink;
                        for jd = ds
							times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
							[ncols, nrows, iv] = myColsAndRowsV(jm, rowstocols);
                        end
                        
                        cclegendstyles = zeros(1,length(times));
                        for jt = 1:length(times)
                            Clines = myLineStyle(length(times), jt, weight);
                            if(isempty(ar.model(jm).plot(jplot).condition))
                                Clines{2} = jeti(jp,:);
                            end
                            
                            countv = 0;
                            for jv = iv
                                countv = countv + 1;
                                [t, v, ~, ~, zero_break] = getDataDoseResponseV(jm, jv, ds, times(jt));
                                if(jp==1)
                                    g = subplot(nrows,ncols,countv);
                                    ar.model(jm).plot(jplot).gv(jv) = g;
                                else
                                    g = ar.model(jm).plot(jplot).gv(jv);
                                end
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, v, Clines{:}, 'LineWidth', linesize);
                                cclegendstyles(ccount) = ltmp;
                                ar.model(jm).data(jd).plot.v(jv,jt) = ltmp;
                                hold(g, 'on');
                                if(jp==np && ~isempty(zero_break))
                                    plot(g, [zero_break zero_break], ylim(g), 'k--');
                                end
                            end
                        end
                    end
                    
                    if(jp==1)
                        figcountv = figcountv + 1;
                    end
                end
            end
        end
%     catch exception
%         fprintf('ERROR for parameter set #%i: %s\n', jp, exception.message);
%     end
end

figcount = 1;
figcountx = 1;
figcountv = 1;
for jm = 1:length(ar.model)
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotYs(jplot) && ar.model(jm).plot(jplot).ny>0)
            h = myRaiseFigure(jm, jplot, ['Y: ' ar.model(jm).plot(jplot).name], figcount);
            
            ny = length(ar.model(jm).data(ar.model(jm).plot(jplot).dLink(1)).y);
            % axis & titles
            for jy = 1:ny
                jd = ar.model(jm).plot(jplot).dLink(1);
                g = ar.model(jm).plot(jplot).gy(jy);
                
                hold(g, 'off');
                mySubplotStyle(g, labelfontsize, labelfonttype);
                
                if(jy == (nrows-1)*ncols + 1)
                    if(~ar.model(jm).plot(jplot).doseresponse)
                        xlabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).tUnits{3}, ar.model(jm).data(jd).tUnits{2}));
                    else
                        xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                    end
                end
                if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
                    ylabel(g, sprintf('log_{10}(%s) [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                else
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                end
                
                if(jy == 1 &&  (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                    if(~ar.model(jm).plot(jplot).doseresponse)
                        legend(g, cclegendstyles, myNameTrafo(ar.model(jm).plot(jplot).condition))
                    else
                        legendtmp = {};
                        ccount = 1;
                        for jt=1:length(times)
                            if(~isempty(conditions))
                                for jc = 1:length(conditions)
                                    legendtmp{ccount} = sprintf('t=%g : %s', times(jt), conditions{jc}); %#ok<AGROW>
                                    ccount = ccount + 1;
                                end
                            else
                                legendtmp{ccount} = sprintf('t=%g', times(jt)); %#ok<AGROW>
                                ccount = ccount + 1;
                            end
                        end
                        legend(g, cclegendstyles, myNameTrafo(legendtmp))
                    end
                end
                
                title(g, myNameTrafo(ar.model(jm).data(jd).y{jy}));
                spacedAxisLimits(g, overplot);
            end
            
            if(saveToFile)
                ar.model(jm).plot(jplot).savePath_FigYMulti = mySaveFigure(h, ar.model(jm).plot(jplot).name, filenameAddition);
            end
            
            figcount = figcount + 1;
        end
        
        if(ar.model(jm).qPlotXs(jplot))
            h = myRaiseFigureX(jm, jplot, ['X: ' ar.model(jm).plot(jplot).name], figcountx);
            [ncols, nrows, nu, nx, iu, ix] = myColsAndRowsX(jm, rowstocols);
            
            % axis & titles
            for jd = ar.model(jm).plot(jplot).dLink
                countu = 0;
                for ju = iu
                    countu = countu + 1;
                    g = ar.model(jm).plot(jplot).gu(ju);
                    
                    hold(g, 'off');
                    
                    title(g, myNameTrafo(ar.model(jm).u{ju}));
                    if(ju == 1 &&  (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            legend(g, cclegendstyles, myNameTrafo(ar.model(jm).plot(jplot).condition))
                        else
                            legendtmp = {};
                            ccount = 1;
                            for jt=1:length(times)
                                if(~isempty(conditions))
                                    for jc = 1:length(conditions)
                                        legendtmp{ccount} = sprintf('t=%g : %s', times(jt), conditions{jc}); %#ok<AGROW>
                                        ccount = ccount + 1;
                                    end
                                else
                                    legendtmp{ccount} = sprintf('t=%g', times(jt)); %#ok<AGROW>
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

                    spacedAxisLimits(g, overplot);
                end
                countx = 0;
                for jx = ix
                    countx = countx + 1;
                    g = ar.model(jm).plot(jplot).gx(jx);
                    hold(g, 'off');
                    
                    title(g, myNameTrafo(ar.model(jm).x{jx}));
                    if(countx+nu == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                        else
                            xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                        end
                    end
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).xUnits{jx,3}, ar.model(jm).xUnits{jx,2}));
                    
                    spacedAxisLimits(g, overplot);
                end
            end
            
            if(saveToFile)
                ar.model(jm).plot(jplot).savePath_FigXMulti = mySaveFigureX(h, ar.model(jm).plot(jplot).name, filenameAddition);
            end
            
            figcountx = figcountx + 1;
        end
        
        if(ar.model(jm).qPlotVs(jplot))
            h = myRaiseFigureV(jm, jplot, ['V: ' ar.model(jm).plot(jplot).name], figcountv);
            [ncols, nrows, iv] = myColsAndRowsV(jm, rowstocols);
            
            % axis & titles
            for jd = ar.model(jm).plot(jplot).dLink
                countv = 0;
                for jv = iv
                    countv = countv + 1;
                    g = ar.model(jm).plot(jplot).gv(jv);
                    hold(g, 'off');
                    
%                         title(g, sprintf('v_{%i}', jv));
%                         fprintf('v%i: %s\n', jv, ar.model(jm).fv{jv});

                    title(g, sprintf('v_{%i}: %s', jv, myNameTrafo(ar.model(jm).fv{jv})));

                    if(countv == (nrows-1)*ncols + 1)
                        if(~ar.model(jm).plot(jplot).doseresponse)
                            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                        else
                            xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                        end
                    end
                    ylabel(g, sprintf('%s [%s]', ar.model(jm).vUnits{jv,3}, ar.model(jm).vUnits{jv,2}));
                    
                    spacedAxisLimits(g, overplot);
                end
            end
            
            if(saveToFile)
                ar.model(jm).plot(jplot).savePath_FigVMulti = mySaveFigureV(h, ar.model(jm).plot(jplot).name, filenameAddition);
            end
            
            figcountv = figcountv + 1;
        end
    end
end

if(exist('hbar','var'))
	close(hbar)
end
ar.p = pReset;



function [t, y, ystd, tExp, yExp, yExpStd] = getData(jm, jd, jy)
global ar
t = ar.model(jm).data(jd).tFine;
y = ar.model(jm).data(jd).yFineSimu(:,jy);
ystd = ar.model(jm).data(jd).ystdFineSimu(:,jy);
if(isfield(ar.model(jm).data(jd), 'yExp'))
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


function [t, u, x, jc] = getDataX(jm, jd)
global ar

if(jd~=0)
    jc = ar.model(jm).data(jd).cLink;
    t = ar.model(jm).data(jd).tFine;
    u = ar.model(jm).condition(jc).uFineSimu(ar.model(jm).data(jd).tLinkFine,:);
    x = ar.model(jm).condition(jc).xFineSimu(ar.model(jm).data(jd).tLinkFine,:);
else
    jc = 1;
    t = ar.model(jm).condition(jc).tFine;
    u = ar.model(jm).condition(jc).uFineSimu;
    x = ar.model(jm).condition(jc).xFineSimu;
end

function [t, v, jc] = getDataV(jm, jd)
global ar

if(jd~=0)
    jc = ar.model(jm).data(jd).cLink;
    t = ar.model(jm).data(jd).tFine;
    v = ar.model(jm).condition(jc).vFineSimu(ar.model(jm).data(jd).tLinkFine,:);
else
    jc = 1;
    t = ar.model(jm).condition(jc).tFine;
    v = ar.model(jm).condition(jc).vFineSimu;
end

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

function [t, v, lb, ub, zero_break] = getDataDoseResponseV(jm, jv, ds, ttime)
global ar

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
            t(ccount,1) = doses(1) - (doses(2)-doses(1)); %#ok<AGROW>
            zero_break = (t(ccount,1)+doses(1))/2;
        end
        v(ccount,1) = ar.model(jm).condition(jc).vExpSimu(jt,jv); %#ok<AGROW>
        
        if(isfield(ar.model(jm).condition(jc), 'vExpUB'))
            lb(ccount,1) = ar.model(jm).condition(jc).vExpLB(jt,jv); %#ok<AGROW>
            ub(ccount,1) = ar.model(jm).condition(jc).vExpUB(jt,jv); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end

%% sub-functions

function C = myLineStyle(n, j, weight)
farben = lines(n);
farben(1,:) = [0 0 0];
C = cell(1,2);
C{1} = 'Color';
C{2} = farben(j,:)*weight + (1-weight);


function h = myRaiseFigure(m, jplot, figname, figcount)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;

if(isfield(ar.model(m).plot(jplot), 'fighandelmulti_y') && ~isempty(ar.model(m).plot(jplot).fighandelmulti_y) && ...
        ar.model(m).plot(jplot).fighandelmulti_y ~= 0 && ...
        sum(ar.model(m).plot(jplot).fighandelmulti_y==openfigs)>0 && ...
		strcmp(get(ar.model(m).plot(jplot).fighandelmulti_y, 'Name'), figname))

	h = ar.model(m).plot(jplot).fighandelmulti_y;
	figure(h);
else
	h = figure('Name', figname, 'NumberTitle','off', ...
		'Units', 'normalized', 'Position', ...
		[0.1+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
	set(h,'Color', figcolor);
	ar.model(m).plot(jplot).fighandelmulti_y = h;
end


function h = myRaiseFigureX(m, jplot, figname, figcount)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;

if(isfield(ar.model(m).plot(jplot), 'fighandelmulti_x') && ~isempty(ar.model(m).plot(jplot).fighandelmulti_x) && ...
        ar.model(m).plot(jplot).fighandelmulti_x ~= 0 && ...
        sum(ar.model(m).plot(jplot).fighandelmulti_x==openfigs)>0 && ...
        strcmp(get(ar.model(m).plot(jplot).fighandelmulti_x, 'Name'), figname))
    
    h = ar.model(m).plot(jplot).fighandelmulti_x;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.4+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.model(m).plot(jplot).fighandelmulti_x = h;
end



function h = myRaiseFigureV(m, jplot, figname, figcount)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot(jplot).time = now;

if(isfield(ar.model(m).plot(jplot), 'fighandelmulti_v') && ~isempty(ar.model(m).plot(jplot).fighandelmulti_v) && ...
        ar.model(m).plot(jplot).fighandelmulti_v ~= 0 && ...
        sum(ar.model(m).plot(jplot).fighandelmulti_v==openfigs)>0 && ...
        strcmp(get(ar.model(m).plot(jplot).fighandelmulti_v, 'Name'), figname))
    
    h = ar.model(m).plot(jplot).fighandelmulti_v;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.8+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.model(m).plot(jplot).fighandelmulti_v = h;
end



function savePath = mySaveFigure(h, name, filenameAddition)
savePath = [arSave '/Figures/Ys' filenameAddition];

if(~exist(savePath, 'dir'))
	mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc', savePath);
eval(['!ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);



function savePath = mySaveFigureX(h, name, filenameAddition)
savePath = [arSave '/Figures/Xs' filenameAddition];

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc2', savePath);
eval(['!ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);


function savePath = mySaveFigureV(h, name, filenameAddition)
savePath = [arSave '/FiguresCI/Vs' filenameAddition];

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



function [ncols, nrows, ny] = myColsAndRows(jm, jd, rowstocols)
global ar
ny = size(ar.model(jm).data(jd).y, 2);
[nrows, ncols] = NtoColsAndRows(ny, rowstocols);


function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);


function [ncols, nrows, nu, nx, iu, ix] = myColsAndRowsX(jm, rowstocols)
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


function [ncols, nrows, iv] = myColsAndRowsV(jm, rowstocols)
global ar
if(~isfield(ar.model(jm), 'qPlotV'))
    nv = size(ar.model(jm).fv, 2);
    iv = 1:nv;
else
    nv = sum(ar.model(jm).qPlotV);
    iv = find(ar.model(jm).qPlotV);
end
[nrows, ncols] = NtoColsAndRows(nv, rowstocols);


function spacedAxisLimits(g, overplot)
[xmin xmax ymin ymax] = axisLimits(g);
xrange = xmax - xmin;
if(xrange == 0)
	xrange = 1;
end
yrange = ymax - ymin;
if(yrange == 0)
	yrange = 1;
end
xlim(g, [xmin-(xrange*overplot) xmax+(xrange*overplot)]);
ylim(g, [ymin-(yrange*overplot) ymax+(yrange*overplot)]);



function [xmin xmax ymin ymax] = axisLimits(g)
p = get(g,'Children');
xmin = nan;
xmax = nan;
ymin = nan;
ymax = nan;
for j = 1:length(p)
	if(~strcmp(get(p(j), 'Type'), 'text'))
		xmin = min([xmin toRowVector(get(p(j), 'XData'))]);
		xmax = max([xmax toRowVector(get(p(j), 'XData'))]);
		%         get(p(j), 'UData')
		%         get(p(j), 'LData')
		%         set(p(j), 'LData', get(p(j), 'LData')*2)
		if(strcmp(get(p(j), 'Type'),'hggroup'))
			ymin = min([ymin toRowVector(get(p(j), 'YData'))-toRowVector(get(p(j), 'LData'))]);
			ymax = max([ymax toRowVector(get(p(j), 'YData'))+toRowVector(get(p(j), 'UData'))]);
		else
			ymin = min([ymin toRowVector(get(p(j), 'YData'))]);
			ymax = max([ymax toRowVector(get(p(j), 'YData'))]);
		end
	end
end



function b = toRowVector(a)
b = a(:)';



