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

matVer = ver('MATLAB');

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
        qDR = ar.model(jm).plot(jplot).doseresponse;
        
        % log 10 dose response axis
        if(isfield(ar.model(jm).plot(jplot), 'doseresponselog10xaxis'))
            logplotting_xaxis = ar.model(jm).plot(jplot).doseresponselog10xaxis;
        else
            logplotting_xaxis = true;
        end
        
        % chi^2, ndata and dr_times
        chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
        ndata = zeros(1,ar.model(jm).plot(jplot).ny);
        dr_times = [];
        for jd = ar.model(jm).plot(jplot).dLink
            if(qDR)
                dr_times = union(dr_times, ar.model(jm).data(jd).tExp); %R2013a compatible
            end
            
            ny = length(ar.model(jm).data(jd).y);
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
        
        ar.model(jm).plot(jplot).chi2 = sum(chi2);
        ar.model(jm).plot(jplot).ndata = sum(ndata);
        
        ar.model(jm).chi2 = ar.model(jm).chi2 + sum(chi2);
        ar.model(jm).ndata = ar.model(jm).ndata + sum(ndata);
        
        if(isempty(dr_times)) % for non dose response
            dr_times = 0;
        end
        
        % conditions
        if(str2double(matVer.Version)>=8.1)
            [conditions, iconditions, jconditions] = ...
                unique(ar.model(jm).plot(jplot).condition,'legacy'); %#ok<ASGLU>
        else
            [conditions, iconditions, jconditions] = ...
                unique(ar.model(jm).plot(jplot).condition); %#ok<ASGLU>
        end
        
        % legends handles and labels
        Clegend = zeros(1,length(dr_times)*length(conditions));
        Clegendlabel = cell(1,length(dr_times)*length(conditions));
        
        if(ar.model(jm).qPlotYs(jplot)==1 && ar.model(jm).plot(jplot).ny>0)

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
            figcount = figcount + 1;
            
            % plotting
            ccount = 1;
            for jt = 1:length(dr_times)
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
                    
                    qUnlog = ar.model(jm).data(jd).logfitting & ...
                        ~ar.model(jm).data(jd).logplotting;
                    qLog = ar.model(jm).data(jd).logplotting;
                    
                    % get data
                    if(qDR)
                        [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, qFit, yExpHl] = ...
                            getDataDoseResponse(jm, ds, dr_times(jt), ...
                            ar.model(jm).plot(jplot).dLink, logplotting_xaxis);
                        y_ssa = [];
                        y_ssa_lb = [];
                        y_ssa_ub = [];
                    else
                        [t, y, ystd, tExp, yExp, yExpStd, lb, ub, yExpHl, ...
                            y_ssa, y_ssa_lb, y_ssa_ub] = getData(jm, jd);
                        qFit = ~isfield(ar.model(jm).data(jd),'qFit') | ...
                            ar.model(jm).data(jd).qFit;
                        zero_break = [];
                    end
                    [tUnits, response_parameter, yLabel, yNames, yUnits] = getInfo(jm, jd);
                    
                    if(isfield(ar.model(jm).data(jd),'plot') && ...
                            isfield(ar.model(jm).data(jd).plot,'y'))
                        hys = ar.model(jm).data(jd).plot.y;
                    else
                        hys = [];
                    end
                    if(isfield(ar.model(jm).data(jd),'plot') && ...
                            isfield(ar.model(jm).data(jd).plot,'y'))
                        hystds = ar.model(jm).data(jd).plot.ystd;
                    else
                        hystds = [];
                    end
                    
                    [hys, hystds, nrows, ncols] = arPlotTrajectories(ccount, ...
                        length(dr_times)*length(jcs), ...
                        t, y, ystd, lb, ub, nfine_dr_plot, ...
                        nfine_dr_method, tExp, yExp, yExpHl, yExpStd, ...
                        y_ssa, y_ssa_lb, y_ssa_ub, ...
                        ar.config.ploterrors, qUnlog, qLog, qFit, ...
                        zero_break, fastPlotTmp, hys, hystds, ...
                        jt==length(dr_times) && jc==jcs(end), qDR, ndata, chi2, ...
                        tUnits, response_parameter, yLabel, yNames, yUnits, ...
                        ar.config.fiterrors, logplotting_xaxis);
                    
                    % line and patch handels
                    ar.model(jm).data(jd).plot.y = hys;
                    ar.model(jm).data(jd).plot.ystd = hystds;
                    
                    % legends
                    Clegend(ccount) = hystds(1);
                    if(qDR)
                        if(~isempty(conditions) && ~isempty(conditions{jc}))
                            Clegendlabel{ccount} = sprintf('t=%g%s : %s', dr_times(jt), ...
                                tUnits{2}, arNameTrafo(conditions{jc}));
                        else
                            Clegendlabel{ccount} = sprintf('t=%g%s', dr_times(jt), ...
                                tUnits{2});
                        end
                    else
                        if(~isempty(conditions) && ~isempty(conditions{jc}))
                            Clegendlabel{ccount} = arNameTrafo(conditions{jc});
                        end
                    end
                    
                    ccount = ccount + 1;
                end
            end
            
            % legend
            if(doLegends && (~isempty(conditions) || qDR))
                g = subplot(nrows, ncols, nrows*ncols);
                axis(g,'off');
                legend(g,Clegend, Clegendlabel);
            end
            
            % optional suptitle
            if(~fastPlotTmp && exist('suptitle','file')==2 && ...
                    isfield(ar.config, 'useSuptitle') && ar.config.useSuptitle)
                suptitle(arNameTrafo([ar.model(jm).name,': ',ar.model(jm).plot(jplot).name]))
            end
            
            % save figure
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
    yExpHl, y_ssa, y_ssa_lb, y_ssa_ub] = getData(jm, jd)

global ar

if(isfield(ar.model(jm).data(jd),'yFineSSA'))
    y_ssa = ar.model(jm).data(jd).yFineSSA;
else
    y_ssa = nan;
end
if(isfield(ar.model(jm).data(jd),'yFineSSA_lb'))
    y_ssa_lb = ar.model(jm).data(jd).yFineSSA_lb;
else
    y_ssa_lb = nan;
end
if(isfield(ar.model(jm).data(jd),'yFineSSA_ub'))
    y_ssa_ub = ar.model(jm).data(jd).yFineSSA_ub;
else
    y_ssa_ub = nan;
end

if(isfield(ar.model(jm).data(jd),'tFine'))
    t = ar.model(jm).data(jd).tFine;
    y = ar.model(jm).data(jd).yFineSimu;
    ystd = ar.model(jm).data(jd).ystdFineSimu;
else
    t = nan;
    y = nan;
    ystd = nan;
end

if(isfield(ar.model(jm).data(jd), 'yExp') && ~isempty(ar.model(jm).data(jd).yExp))
    tExp = ar.model(jm).data(jd).tExp;
    yExp = ar.model(jm).data(jd).yExp;
    if(ar.config.fiterrors == -1)
        yExpStd = ar.model(jm).data(jd).yExpStd;
    else
        if(isfield(ar.model(jm).data(jd),'ystdExpSimu'))
            yExpStd = ar.model(jm).data(jd).ystdExpSimu;
        else
            yExpStd = nan;
        end
    end
    if(isfield(ar.model(jm).data(jd),'highlight'))
        hl = ar.model(jm).data(jd).highlight;
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
    lb = ar.model(jm).data(jd).yFineLB;
    ub = ar.model(jm).data(jd).yFineUB;
else
    lb = [];
    ub = [];
end




function [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, data_qFit, yExpHl] = ...
    getDataDoseResponse(jm, ds, ttime, dLink, logplotting_xaxis)
global ar


zero_break = [];
data_qFit = true;

ccount = 1;
for jd = ds
    data_qFit = data_qFit & ar.model(jm).data(jd).qFit;
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
        y(ccount,:) = ar.model(jm).data(jd).yFineSimu(jtfine,:); %#ok<AGROW>
        ystd(ccount,:) = ar.model(jm).data(jd).ystdFineSimu(jtfine,:); %#ok<AGROW>
        
        yExp(ccount,:) = ar.model(jm).data(jd).yExp(jt,:); %#ok<AGROW>
        yExpHl(ccount,:) = NaN; %#ok<AGROW>
        if(isfield(ar.model(jm).data(jd),'highlight'))
            if(ar.model(jm).data(jd).highlight(jt,:)~=0)
                yExpHl(ccount,:) = yExp(ccount,:);                 %#ok<AGROW>
            end
        end
        if(ar.config.fiterrors == -1)
            yExpStd(ccount,:) = ar.model(jm).data(jd).yExpStd(jt,:); %#ok<AGROW>
        else
            yExpStd(ccount,:) = ar.model(jm).data(jd).ystdExpSimu(jt,:); %#ok<AGROW>
        end
        if(isfield(ar.model(jm).data(jd), 'yExpUB'))
            lb(ccount,:) = ar.model(jm).data(jd).yFineLB(jtfine,:); %#ok<AGROW>
            ub(ccount,:) = ar.model(jm).data(jd).yFineUB(jtfine,:); %#ok<AGROW>
        else
            lb = [];
            ub = [];
        end
        
        ccount = ccount + 1;
    end
end

[tExp,itexp] = sort(tExp);
yExp = yExp(itexp,:);
yExpHl = yExpHl(itexp,:);

[t,it] = sort(t);
y = y(it,:);
ystd = ystd(it,:);
if(~isempty(lb))
    lb = lb(it,:);
    ub = ub(it,:);
end


function [tUnits, response_parameter, yLabel, yNames, yUnits] = getInfo(jm, jd)
global ar
tUnits = ar.model(jm).data(jd).tUnits;
response_parameter = ar.model(jm).data(jd).response_parameter;
yLabel = ar.model(jm).data(jd).y;
if(isfield(ar.model(jm).data(jd), 'yNames'))
    yNames = ar.model(jm).data(jd).yNames;
else
    yNames = [];
end 
yUnits = ar.model(jm).data(jd).yUnits;

