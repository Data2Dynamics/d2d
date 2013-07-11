% Plot models V
%
% arPlotV(saveToFile, fastPlot)
%
% saveToFile    [false]
% fastPlot      [false]

function arPlotV(saveToFile, fastPlot)

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
        if(ar.model(jm).qPlotVs(jplot)==1)
            if(~ar.model(jm).isReactionBased)
                fprintf('arPlotV: model %i is not reaction based, plotting omitted\n', jm);
                return;
            else
                if(ar.config.ploterrors == -1)
                    [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['CI-V: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
                else
                    [h, fastPlotTmp] = myRaiseFigure(jm, jplot, ['V: ' ar.model(jm).plot(jplot).name], figcount, fastPlot);
                end
                
                % plotting
                ccount = 1;
                if(~ar.model(jm).plot(jplot).doseresponse)
                    cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                    
                    for jd = ar.model(jm).plot(jplot).dLink
                        [t, v, vlb, vub, jc] = getData(jm, jd);
                        
                        % rows and cols
                        [ncols, nrows, iv] = myColsAndRows(jm, rowstocols);
                        
                        Clines = myLineStyle(length(ar.model(jm).plot(jplot).dLink), ccount);
                        countv = 0;
                        for jv = iv
                            countv = countv + 1;
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,countv);
                                ar.model(jm).plot(jplot).gv(jv) = g;
                                mySubplotStyle(g, labelfontsize, labelfonttype);
                                ltmp = plot(g, t, v(:,jv), Clines{:});
                                cclegendstyles(ccount) = ltmp;
                                if(jd~=0)
                                    ar.model(jm).data(jd).plot.v(jv,jc) = ltmp;
                                else
                                    ar.model(jm).condition(jc).plot.v(jv,jc) = ltmp;
                                end
                                hold(g, 'on');
                                if(ar.config.ploterrors == -1)
                                    tmpx = [t(:); flipud(t(:))];
                                    tmpy = [vub(:,jv); flipud(vlb(:,jv))];
                                    ltmp = patch(tmpx, tmpy, tmpx*0-2*eps, 'r');
                                    set(ltmp, 'FaceColor', Clines{2}*0.1+0.9, 'EdgeColor', Clines{2}*0.1+0.9);
                                    ltmp2 = patch(tmpx, tmpy, tmpx*0-eps, 'r');
                                    set(ltmp2, 'LineStyle', '--', 'FaceColor', 'none', 'EdgeColor', Clines{2}*0.3+0.7);
                                end
                            else
                                if(jd~=0)
                                    set(ar.model(jm).data(jd).plot.v(jv,jc), 'YData', v(:,jv));
                                else
                                    set(ar.model(jm).condition(jc).plot.v(jv,jc), 'YData', v(:,jv));
                                end
                            end
                        end
                        ccount = ccount + 1;
                    end
                else
                    times = [];
                    for jd = ar.model(jm).plot(jplot).dLink
						times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                        [ncols, nrows, iv] = myColsAndRows(jm, rowstocols);
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
                            
                            countv = 0;
                            for jv = iv
                                countv = countv + 1;
                                [t, v, lb, ub, zero_break] = getDataDoseResponseV(jm, jv, ds, times(jt));
                                if(~fastPlotTmp)
                                    g = subplot(nrows,ncols,countv);
                                    ar.model(jm).plot(jplot).gv(jv) = g;
                                    mySubplotStyle(g, labelfontsize, labelfonttype);
                                    ltmp = plot(g, t, v, Clines{:});
                                    cclegendstyles(ccount) = ltmp;
                                    ar.model(jm).data(jd).plot.v(jv,jt,jc) = ltmp;
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
                                    set(ar.model(jm).data(jd).plot.v(jv,jt,jc), 'YData', v);
                                end
                            end
                            ccount = ccount + 1;
                        end
                    end
                end
                
                % axis & titles
                for jd = ar.model(jm).plot(jplot).dLink
                    countv = 0;
                    for jv = iv
                        countv = countv + 1;
                        g = ar.model(jm).plot(jplot).gv(jv);
                        if(~fastPlotTmp)
                            hold(g, 'off');
                            
                            title(g, sprintf('v_{%i}', jv));
                            %                         fprintf('v%i: %s\n', jv, ar.model(jm).fv{jv});
                            %                         title(g, sprintf('v_{%i}: %s', jv, myNameTrafo(ar.model(jm).fv{jv})));
                            
                            if(countv == (nrows-1)*ncols + 1)
                                if(~ar.model(jm).plot(jplot).doseresponse)
                                    xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                                else
                                    xlabel(g, sprintf('log_{10}(%s)', myNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                                end
                            end
                            ylabel(g, sprintf('%s [%s]', ar.model(jm).vUnits{jv,3}, ar.model(jm).vUnits{jv,2}));
                        end
                        arSpacedAxisLimits(g, overplot);
                    end
                end
                
                if(saveToFile)
                    if(ar.config.ploterrors == -1)
                        ar.model(jm).plot(jplot).savePath_FigVCI = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                    else
                        ar.model(jm).plot(jplot).savePath_FigV = mySaveFigure(h, ar.model(jm).plot(jplot).name);
                    end
                end
                
                figcount = figcount + 1;
            end
        else
            try %#ok<TRYNC>
                close(ar.model(jm).plot(jplot).fighandel_v)
            end
            ar.model(jm).plot(jplot).fighandel_v = [];
        end
    end
end

function [t, v, vlb, vub, jc] = getData(jm, jd)
global ar

if(jd~=0)
    jc = ar.model(jm).data(jd).cLink;
    t = ar.model(jm).data(jd).tFine;
    v = ar.model(jm).condition(jc).vFineSimu(ar.model(jm).data(jd).tLinkFine,:);
    if(isfield(ar.model(jm).condition(jc), 'vExpUB'))
        vlb = ar.model(jm).condition(jc).vFineLB(ar.model(jm).data(jd).tLinkFine,:);
        vub = ar.model(jm).condition(jc).vFineUB(ar.model(jm).data(jd).tLinkFine,:);
    else
        vlb = [];
        vub = [];
    end
else
    jc = 1;
    t = ar.model(jm).condition(jc).tFine;
    v = ar.model(jm).condition(jc).vFineSimu;
    if(isfield(ar.model(jm).condition(jc), 'vExpUB'))
        vlb = ar.model(jm).condition(jc).vFineLB;
        vub = ar.model(jm).condition(jc).vFineUB;
    else
        vlb = [];
        vub = [];
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
    if(isfield(ar.model(m).plot(jplot), 'fighandel_vCI') && ~isempty(ar.model(m).plot(jplot).fighandel_vCI) && ...
            ar.model(m).plot(jplot).fighandel_vCI ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_vCI==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_vCI, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_vCI;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.8+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_vCI = h;
        fastPlotTmp = false;
    end    
else
    if(isfield(ar.model(m).plot(jplot), 'fighandel_v') && ~isempty(ar.model(m).plot(jplot).fighandel_v) && ...
            ar.model(m).plot(jplot).fighandel_v ~= 0 && ...
            sum(ar.model(m).plot(jplot).fighandel_v==openfigs)>0 && ...
            strcmp(get(ar.model(m).plot(jplot).fighandel_v, 'Name'), figname))
        
        h = ar.model(m).plot(jplot).fighandel_v;
        if(~fastPlot)
            figure(h);
        end
    else
        h = figure('Name', figname, 'NumberTitle','off', ...
            'Units', 'normalized', 'Position', ...
            [0.8+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
        set(h,'Color', figcolor);
        ar.model(m).plot(jplot).fighandel_v = h;
        fastPlotTmp = false;
    end
end
if(~fastPlot)
    clf
end



function savePath = mySaveFigure(h, name)
global ar
if(ar.config.ploterrors == -1)
    savePath = [arSave '/FiguresCI/V'];
else
    savePath = [arSave '/Figures/V'];
end

if(~exist(savePath, 'dir'))
    mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc2', savePath);
eval(['!ps2pdf  -dEPSCrop -dAutoRotatePages=/None ' savePath '.eps '  savePath '.pdf']);



function str = mypath(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');



function str = myNameTrafo(str)
str = strrep(str, '_', '\_');



function mySubplotStyle(g, labelfontsize, labelfonttype)
set(g, 'FontSize', labelfontsize);
set(g, 'FontName', labelfonttype);



function [ncols, nrows, iv] = myColsAndRows(jm, rowstocols)
global ar
if(~isfield(ar.model(jm), 'qPlotV'))
    nv = size(ar.model(jm).fv, 2);
    iv = 1:nv;
else
    nv = sum(ar.model(jm).qPlotV);
    iv = find(ar.model(jm).qPlotV);
end
[nrows, ncols] = NtoColsAndRows(nv, rowstocols);



function [nrows, ncols] = NtoColsAndRows(n, rowstocols)
nrows = ceil(n^rowstocols);
ncols = ceil(n / nrows);



