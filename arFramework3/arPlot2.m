% Plot models and datasets
%
% arPlot(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics)
%
% saveToFile    [false]
% fastPlot      [false]
% silent        [false]
% evalfun       [true]
% doLegends     [true]
% dynamics:     [true]

function arPlot2(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics)

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
if(~exist('silent','var'))
    silent = false;
end
if(~exist('evalfun','var'))
    evalfun = true;
end
if(~exist('doLegends','var'))
    doLegends = true;
end
if(~exist('dynamics','var'))
    dynamics = true;
end

if(evalfun)
    try
        arSimu(false, true, dynamics);
        if(silent)
            arChi2(false, [], dynamics);
        else
            arChi2;
        end
    catch err_id
        if(silent)
            disp(err_id.message);
        end
    end
end

if(~isfield(ar.model, 'qPlotYs'))
    for jm=1:length(ar.model)
        if(length(ar.model(jm).plot) > 5)
            fprintf('Automatic plotting disabled for model %i. Please use arTuner for plotting.\n', jm);
            ar.model(jm).qPlotYs = false(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotXs = false(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotVs = false(1,length(ar.model(jm).plot));
        else
            ar.model(jm).qPlotYs = true(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotXs = false(1,length(ar.model(jm).plot));
            ar.model(jm).qPlotVs = false(1,length(ar.model(jm).plot));
        end
    end
end

matVer = ver('MATLAB');

if(~exist('saveToFile','var'))
    saveToFile = false;
end
if(~exist('fastPlot','var'))
    fastPlot = false;
end
if(~exist('doLegends','var'))
    doLegends = true;
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
    ar.model(jm).chi2 = 0;
    ar.model(jm).ndata = 0;
    
    for jplot = 1:length(ar.model(jm).plot)
        if(ar.model(jm).qPlotYs(jplot)==1 && ar.model(jm).plot(jplot).ny>0)

            qDR = ar.model(jm).plot(jplot).doseresponse;
            
            % setup figure
            if(ar.config.ploterrors == -1)
                [h, fastPlotTmp] = arRaiseFigure(ar.model(jm).plot(jplot), ...
                    'fighandel_yCI', ['CI-Y: ' ar.model(jm).plot(jplot).name], ...
                    figcount, fastPlot);
                ar.model(jm).plot(jplot).fighandel_yCI = h;
            else
                [h, fastPlotTmp] = arRaiseFigure(ar.model(jm).plot(jplot), ...
                    'fighandel_y', ['Y: ' ar.model(jm).plot(jplot).name], ...
                    figcount, fastPlot);
                ar.model(jm).plot(jplot).fighandel_y = h;
            end
            
            % log 10 dose response axis
            if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
                logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
            else
                logplotting_xaxis = true;
            end
            
            % chi^2, ndata and times
            chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
            ndata = zeros(1,ar.model(jm).plot(jplot).ny);
            times = [];
            for jd = ar.model(jm).plot(jplot).dLink
                if(qDR)
                    times = union(times, ar.model(jm).data(jd).tExp); %R2013a compatible
                end
                
                ny = size(ar.model(jm).data(jd).y, 2);
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
            if(isempty(times)) % for non dose response
                times = 0;
            end
            
            % conditions
            if(str2double(matVer.Version)>=8.1)
                [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition,'legacy'); %#ok<ASGLU>
            else
                [conditions, iconditions, jconditions] = unique(ar.model(jm).plot(jplot).condition); %#ok<ASGLU>
            end
            
            % legends handles
            cclegendstyles = zeros(1,length(times)*length(conditions));
            
            % plotting
            ccount = 1;
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
                    
                    % rows and cols
                    ny = size(ar.model(jm).data(jd).y, 2);
                    [nrows, ncols] = arNtoColsAndRows(ny);
                    if(nrows*ncols == ny)
                        [nrows, ncols] = arNtoColsAndRows(ny+1);
                    end
                    
                    % styles
                    Clines = arLineMarkersAndColors(ccount, ...
                        length(times)*length(jcs), ...
                        [], 'none', '-');
                    ClinesExp = arLineMarkersAndColors(ccount, ...
                        length(times)*length(jcs), ...
                        [], 'none', 'none');
                    
                    for jy = 1:ny
                        if(qDR)
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, qFit, yExpHl] = ...
                                getDataDoseResponse(jm, jy, ds, times(jt), ...
                                ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                            y_ssa = [];
                            y_ssa_lb = [];
                            y_ssa_ub = [];
                        else
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, yExpHl, ...
                                y_ssa, y_ssa_lb, y_ssa_ub] = getData(jm, jd, jy);
                            qFit = ~isfield(ar.model(jm).data(jd),'qFit') || ...
                                ar.model(jm).data(jd).qFit(jy);
                            zero_break = [];
                        end
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
                        
                        qUnlog = ar.model(jm).data(jd).logfitting(jy) && ...
                            ~ar.model(jm).data(jd).logplotting(jy);
                        
                        if(~fastPlotTmp)
                            g = subplot(nrows,ncols,jy);
                            hold(g, 'on');
                            ar.model(jm).plot(jplot).gy(jy) = g;
                            
                            % call arPlotTrajectory
                            [hy, hystd] = arPlotTrajectory(t, y, ystd, lb, ub, ...
                                tExp, yExp, yExpHl, yExpStd, ...
                                y_ssa, y_ssa_lb, y_ssa_ub, ...
                                ar.config.ploterrors, Clines, ClinesExp, qUnlog, ...
                                [], [], qFit, zero_break);
                            
                            % save handles for fast plotting
                            ar.model(jm).data(jd).plot.y(jy) = hy;
                            if(isempty(hystd))
                                hystd = hy;
                            end
                            cclegendstyles(ccount) = hystd;
                            ar.model(jm).data(jd).plot.ystd(jy) = hystd;
                        else
                            % call arPlotTrajectory
                            arPlotTrajectory(t, y, ystd, lb, ub, ...
                                tExp, yExp, yExpHl, yExpStd, ...
                                y_ssa, y_ssa_lb, y_ssa_ub, ...
                                ar.config.ploterrors, Clines, ClinesExp, qUnlog, ...
                                ar.model(jm).data(jd).plot.y(jy), ...
                                ar.model(jm).data(jd).plot.ystd(jy), qFit, []);
                        end
                    end
                    ccount = ccount + 1;
                end
            end
            
            % axis & titles
            
            % optional suptitle
            if(~fastPlotTmp && exist('suptitle','file')==2 && ...
                    isfield(ar.config, 'useSuptitle') && ar.config.useSuptitle)
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
                    if(ar.model(jm).data(jd).logfitting(jy) && ar.model(jm).data(jd).logplotting(jy))
                        ylabel(g, sprintf('log_{10}(%s) [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    else
                        ylabel(g, sprintf('%s [%s]', ar.model(jm).data(jd).yUnits{jy,3}, ar.model(jm).data(jd).yUnits{jy,2}));
                    end
                    
                    if(doLegends && jy == 1 && (~isempty(ar.model(jm).plot(jplot).condition) || ar.model(jm).plot(jplot).doseresponse))
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
                            legloc(1:2) = grefloc(1:2) + [grefloc(3)-legloc(3) 0];
                            % legloc(1:2) = [1-legloc(3) 0];
                            set(hl, 'Position', legloc);
                        end
                    end
                end
                titstr = {};
                if(isfield(ar.model(jm).data(jd), 'yNames') && ~isempty(ar.model(jm).data(jd).yNames{jy}) && ...
                        ~strcmp(ar.model(jm).data(jd).yNames{jy}, ar.model(jm).data(jd).y{jy}))
                    titstr{1} = [arNameTrafo(ar.model(jm).data(jd).yNames{jy}) ' (' arNameTrafo(ar.model(jm).data(jd).y{jy}) ')'];
                else
                    titstr{1} = arNameTrafo(ar.model(jm).data(jd).y{jy});
                end
                if(isfield(ar.model(jm).data(jd), 'yExp'))
                    if(ndata(jy)>0)
                        if(ar.config.fiterrors == 1)
                            titstr{2} = sprintf('-2 log(L)_{%i} = %g', ndata(jy), 2*ndata(jy)*log(sqrt(2*pi)) + chi2(jy));
                        else
                            titstr{2} = sprintf('chi^2_{%i} = %g', ndata(jy), chi2(jy));
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
            
            if(saveToFile)
                if(ar.config.ploterrors == -1)
                    ar.model(jm).plot(jplot).savePath_FigYCI = arSaveFigure(h, ...
                        ar.model(jm).plot(jplot).name, '/FiguresCI/Y');
                else
                    ar.model(jm).plot(jplot).savePath_FigY = arSaveFigure(h, ...
                        ar.model(jm).plot(jplot).name, '/Figures/Y');
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



function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, ...
    yExpHl, y_ssa, y_ssa_lb, y_ssa_ub] = getData(jm, jd, jy)

global ar

if(isfield(ar.model(jm).data(jd),'yFineSSA'))
    y_ssa = ar.model(jm).data(jd).yFineSSA(:,jy,:);
else
    y_ssa = nan;
end
if(isfield(ar.model(jm).data(jd),'yFineSSA_lb'))
    y_ssa_lb = ar.model(jm).data(jd).yFineSSA_lb(:,jy,:);
else
    y_ssa_lb = nan;
end
if(isfield(ar.model(jm).data(jd),'yFineSSA_ub'))
    y_ssa_ub = ar.model(jm).data(jd).yFineSSA_ub(:,jy,:);
else
    y_ssa_ub = nan;
end

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
    if(ar.config.fiterrors == -1)
        yExpStd = ar.model(jm).data(jd).yExpStd(:,jy);
    else
        if(isfield(ar.model(jm).data(jd),'ystdExpSimu'))
            yExpStd = ar.model(jm).data(jd).ystdExpSimu(:,jy);
        else
            yExpStd = nan;
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
        yExpHl(ccount,1) = NaN; %#ok<AGROW>
        if(isfield(ar.model(jm).data(jd),'highlight'))
            if(ar.model(jm).data(jd).highlight(jt,jy)~=0)
                yExpHl(ccount,1) = yExp(ccount,1);                 %#ok<AGROW>
            end
        end
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

[tExp,itexp] = sort(tExp);
yExp = yExp(itexp);
yExpHl = yExpHl(itexp);

[t,it] = sort(t);
y = y(it);
ystd = ystd(it);
if(~isempty(lb))
    lb = lb(it);
    ub = ub(it);
end


