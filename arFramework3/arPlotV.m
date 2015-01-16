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
if(ar.config.ploterrors == -1)
    linesize = 0.5;
else
    linesize = 2;
end

if(isfield(ar.config,'nfine_dr_plot'))
    nfine_dr_plot = ar.config.nfine_dr_plot;
    nfine_dr_method = ar.config.nfine_dr_method;
else
    nfine_dr_plot = 1;
    nfine_dr_method = 'spline';
end

figcount = 1;
for jm = 1:length(ar.model)
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotVs(jplot)==1)
            if(~ar.model(jm).isReactionBased)
                fprintf('arPlotV: model %i is not reaction based, plotting omitted\n', jm);
                return;
            else
                if(ar.config.ploterrors == -1)
                    [h, fastPlotTmp] = arRaiseFigure(ar.model(jm).plot(jplot), ...
                        'fighandel_vCI', ['CI-V: ' ar.model(jm).plot(jplot).name], ...
                        figcount, fastPlot, 3);
                else
                    [h, fastPlotTmp] = arRaiseFigure(ar.model(jm).plot(jplot), ...
                        'fighandel_v', ['V: ' ar.model(jm).plot(jplot).name], ...
                        figcount, fastPlot, 3);
                end
                
                if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
                    logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
                else
                    logplotting_xaxis = true;
                end
                
                % plotting
                ccount = 1;
                if(~ar.model(jm).plot(jplot).doseresponse)
                    cclegendstyles = zeros(1,length(ar.model(jm).plot(jplot).dLink));
                    
                    for jd = ar.model(jm).plot(jplot).dLink
                        [t, v, vlb, vub, jc] = getData(jm, jd);
                        
                        % rows and cols
                        [ncols, nrows, iv] = myColsAndRows(jm);
                        
                        Clines = arLineMarkersAndColors(ccount, ...
                            length(ar.model(jm).plot(jplot).dLink), ...
                            [], 'none', '-');
                        
                        countv = 0;
                        for jv = iv
                            countv = countv + 1;
                            if(~fastPlotTmp)
                                g = subplot(nrows,ncols,countv);
                                ar.model(jm).plot(jplot).gv(jv) = g;
                                arSubplotStyle(g);
                                ltmp = plot(g, t, v(:,jv), Clines{:}, 'LineWidth', linesize);
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
                        [ncols, nrows, iv] = myColsAndRows(jm);
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
                            
                            countv = 0;
                            for jv = iv
                                countv = countv + 1;
                                [t, v, lb, ub, zero_break] = getDataDoseResponseV(jm, jv, ds, times(jt), ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                                if(length(unique(t))==1)
                                    t = [t-0.1; t+0.1];
                                    v = [v; v]; %#ok<AGROW>
                                    lb = [lb; lb]; %#ok<AGROW>
                                    ub = [ub; ub]; %#ok<AGROW>
                                elseif(nfine_dr_plot>10)
                                    tf = linspace(min(t), max(t), nfine_dr_plot);
                                    v = interp1(t,v,tf,nfine_dr_method);
                                    if(~isempty(lb))
                                        lb = interp1(t,lb,tf,nfine_dr_method);
                                        ub = interp1(t,ub,tf,nfine_dr_method);
                                    end
                                    t = tf;
                                end
                                
                                if(~fastPlotTmp)
                                    g = subplot(nrows,ncols,countv);
                                    ar.model(jm).plot(jplot).gv(jv) = g;
                                    arSubplotStyle(g);
                                    ltmp = plot(g, t, v, Clines{:}, 'LineWidth', linesize);
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
                            if(isempty(ar.model(jm).v{jv}))
                                title(g, sprintf('v_{%i}', jv));
                            else
                                title(g,arNameTrafo(ar.model(jm).v{jv}));
                            end
                            fprintf('v%i: %s\n', jv, ar.model(jm).fv{jv});
                            % title(g, sprintf('v_{%i}: %s', jv, arNameTrafo(ar.model(jm).fv{jv})));
                            
                            if(countv == (nrows-1)*ncols + 1)
                                if(~ar.model(jm).plot(jplot).doseresponse)
                                    xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
                                else
                                    xlabel(g, sprintf('log_{10}(%s)', arNameTrafo(ar.model(jm).data(jd).condition(1).parameter)));
                                end
                            end
                            % ylabel(g, sprintf('%s [%s]', ar.model(jm).vUnits{jv,3}, ar.model(jm).vUnits{jv,2}));
                        end
                        arSpacedAxisLimits(g);
                    end
                end
                
                if(saveToFile)
                    if(ar.config.ploterrors == -1)
                        ar.model(jm).plot(jplot).savePath_FigVCI = arSaveFigure(h, ...
                            ar.model(jm).plot(jplot).name, '/FiguresCI/V');
                    else
                        ar.model(jm).plot(jplot).savePath_FigV = arSaveFigure(h, ...
                            ar.model(jm).plot(jplot).name, '/Figures/V');
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



function [t, v, lb, ub, zero_break] = getDataDoseResponseV(jm, jv, ds, ttime, dLink, logplotting_xaxis)
global ar

zero_break = [];

ccount = 1;
for jd = ds
    for jc = 1:length(ar.model(jm).data(jd).condition)
        if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, ar.model(jm).data(jd).response_parameter))
            jcondi = jc;
        end
    end
    
    jc = ar.model(jm).data(jd).cLink;
    qt = ar.model(jm).condition(jc).tExp == ttime;
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
[t,it] = sort(t);
v = v(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end


%% sub-functions



function [ncols, nrows, iv] = myColsAndRows(jm)
global ar
if(~isfield(ar.model(jm), 'qPlotV'))
    nv = size(ar.model(jm).fv, 2);
    iv = 1:nv;
else
    nv = sum(ar.model(jm).qPlotV);
    iv = find(ar.model(jm).qPlotV);
end
[nrows, ncols] = arNtoColsAndRows(nv);




