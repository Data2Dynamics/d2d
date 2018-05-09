% Plot models and datasets (new function)
%
% hs = arPlot2(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics, hs)
%
% hs: figure handles;
%
% saveToFile    [false]
% fastPlot      [false]
% silent        [false]
% evalfun       [true]
% doLegends     [true]
% dynamics:     [true]
% hs:           []      custom figure handels
% 
% 
% Doku: 
% https://github.com/Data2Dynamics/d2d/wiki/Plotting-options-and-the-meaning-of-ar.config.ploterrors

function varargout = arPlot2(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics, hs)

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
if(~exist('hs','var'))
	hs = [];
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

matVer = ver('MATLAB');

if(evalfun)
    try
        arSimu(false, true, dynamics);
    catch err_id
        if(~silent)
            disp(err_id.message);
        end
    end
    try
        if(silent)
            arCalcMerit(false, [], dynamics);
        else
            arCalcMerit;
        end
    catch err_id
        if(~silent)
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
hsnew = [];

fiterrors = ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0);

figcount = 1;
for jm = 1:length(ar.model)
    ar.model(jm).chi2 = 0;
    ar.model(jm).ndata = 0;
    
    for jplot = 1:length(ar.model(jm).plot)
        qDR = ar.model(jm).plot(jplot).doseresponse;

        % Determine transformation of the independent axis
        [xtrafo, xLabel] = arGetPlotXTrafo(jm, jplot);
        
        % chi^2, ndata and dr_times
        chi2 = zeros(1,ar.model(jm).plot(jplot).ny);
        ndata = zeros(1,ar.model(jm).plot(jplot).ny);
        dr_times = [];
        for jd = ar.model(jm).plot(jplot).dLink
            if(isfield(ar.model(jm),'data'))
                if(qDR)
                    dr_times = union(dr_times, ar.model(jm).data(jd).tExp); %R2013a compatible
                end
                
                ny = length(ar.model(jm).data(jd).y);
                for jy = 1:ny
                    % chi^2 & ndata
                    if(ar.model(jm).data(jd).qFit(jy)==1)
                        chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2(jy);
                        ndata(jy) = ndata(jy) + ar.model(jm).data(jd).ndata(jy);
                        if fiterrors == 1
                            chi2(jy) = chi2(jy) + ar.model(jm).data(jd).chi2err(jy);
                        end
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
            [conditions, ~, jconditions] = ...
                unique(ar.model(jm).plot(jplot).condition,'legacy');
        else
            [conditions, ~, jconditions] = ...
                unique(ar.model(jm).plot(jplot).condition);
        end
        
        fighandel_name = {'fighandel_y', 'fighandel_x', 'fighandel_v'};
        fig_name = {'Y: ', 'X: ', 'V: '};
        linehandle_name = {'y','x','v'};
        savepath_name = {'Y','X','V'};
        qplotname = {'qPlotYs', 'qPlotXs', 'qPlotVs'};
        
        didPlot = false;
        
        % jype = 1: data y is plotted
        % jype = 2: dynamics x is plotted
        % jype = 3: v is plotted
        if(isfield(ar.model(jm),'data'))
            jtypes = 1:3;
        else
            jtypes = 2:3;
        end
        
        for jtype = jtypes
            if jtype ==1
                plotopt_str = ['_ploterrors' num2str(ar.config.ploterrors) '_fiterrors' num2str(ar.config.fiterrors)];
            else
                plotopt_str = '';
            end
            
            % legends handles and labels
            Clegend = zeros(1,length(dr_times)*length(conditions));
            Clegendlabel = cell(1,length(dr_times)*length(conditions));
            
            if(ar.model(jm).(qplotname{jtype})(jplot)==1 && (jtype~=1 || ar.model(jm).plot(jplot).ny>0))
                didPlot = true;
                
                % setup figure
                [h, fastPlotTmp] = arRaiseFigure(ar.model(jm).plot(jplot), ...
                    [fighandel_name{jtype} ], [fig_name{jtype} ar.model(jm).plot(jplot).name plotopt_str], ...
                    figcount, fastPlot, jtype, hs);
                ar.model(jm).plot(jplot).([fighandel_name{jtype}]) = h;
                hsnew(end+1) = h; %#ok<AGROW>
                
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
                        
                        % get data
                        if(qDR)
                            
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, zero_break, qFit, yExpHl] = ...
                                arGetDataDoseResponse(jm, ds, dr_times(jt), ...
                                ar.model(jm).plot(jplot).dLink, jtype, xtrafo);                            
                            
                            plotopt = NaN(1,size(y,2));
                            if jtype ==1
                                for jy=1:size(y,2)
                                    plotopt(jy) = arWhichYplot(jm,ds,[],jy);
                                end
                            else
                                for jy=1:size(y,2)
                                    plotopt(jy) = 1;
                                end
                            end
                            
                            % TODO
                            t_ppl = [];
                            y_ppl_ub = [];
                            y_ppl_lb = [];
                            
                            % TODO
                            y_ssa = [];
                            y_ssa_lb = [];
                            y_ssa_ub = [];
                            dydt = [];
                            t_ppl = [];
                            y_ppl_ub = [];
                            y_ppl_lb = [];
                        else
                            [t, y, ystd, tExp, yExp, yExpStd, lb, ub, yExpHl, dydt, ...
                                y_ssa, y_ssa_lb, y_ssa_ub, qFit, t_ppl, y_ppl_ub, y_ppl_lb] = arGetData(jm, jd, jtype);
                            plotopt = NaN(1,size(y,2));
                            if jtype ==1
                                for jy=1:size(y,2)
                                    plotopt(jy) = arWhichYplot(jm,jd,[],jy);
                                end
                            else
                                for jy=1:size(y,2)
                                    plotopt(jy) = 1;
                                    if(~isempty(t_ppl))
                                        plotopt(jy) = 4;
                                    end
                                end
                            end
                            zero_break = [];
                        end
                        
                        [tUnits, response_parameter, titles, yNames, yLabel, iy, ...
                            hys, hystds, hysss] = ...
                            arGetInfo(jm, jd, jtype, linehandle_name{jtype});

                        % Plotting of observables
                        if(jtype==1)
                            % Central point where the transformations are handled.
                            [trafos, yLabel] = arGetPlotYTrafo(jm, jd, jplot);
                        else
                            trafos = cell(1, numel(iy));
                            for jy = 1 : numel(iy)
                                trafos{jy} = @(x) x;
                            end
                        end
                        
                        % call arPlotTrajectories
                        [hys, hystds, hysss, nrows, ncols] = arPlotTrajectories(ccount, ...
                            length(dr_times)*length(jcs), ...
                            t, y, ystd, lb, ub, nfine_dr_plot, ...
                            nfine_dr_method, tExp, yExp, yExpHl, yExpStd, ...
                            y_ssa, y_ssa_lb, y_ssa_ub, plotopt, trafos, qFit, zero_break, ...
                            fastPlotTmp, hys, hystds, hysss, dydt, ...
                            jt==length(dr_times) && jc==jcs(end), ndata, chi2, ...
                            titles, yNames, xLabel, yLabel, fiterrors, iy, t_ppl, y_ppl_ub, y_ppl_lb, ...
                            ar.config.atol);
        
                        % save handels
                        if(jd~=0)
                            ar.model(jm).data(jd).plot.(linehandle_name{jtype}) = hys;
                            if(jtype == 1)
                                ar.model(jm).data(jd).plot.ystd = hystds;
                            end
                            if(jtype == 2)
                                ar.model(jm).data(jd).plot.xss = hysss;
                            end
                        else
                            ar.model(jm).plot.(linehandle_name{jtype}) = hys;
                            if(jtype == 1)
                                ar.model(jm).plot.ystd = hystds;
                            end
                            if(jtype == 2)
                                ar.model(jm).plot.xss = hysss;
                            end
                        end

                        
                        % legends
                        if(~isempty(hys))
                            if(jtype == 1)
                                inonzero = find(hystds~=0);
                                Clegend(ccount) = hystds(inonzero(1));
                            else
                                inonzero = find(hys~=0);
                                Clegend(ccount) = hys(inonzero(1));
                            end
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
                        end
                        
                        ccount = ccount + 1;
                    end
                end
                
                % legend
                if(doLegends && ~fastPlot && (~isempty(conditions) && nrows*ncols>0 || qDR))
                    % The following 6 commented lines got really slow in
                    % MATLAB R2016a, don't know why, should be revised at
                    % some point
%                     g = subplot(nrows, ncols, nrows*ncols);
%                     lpos = get(g,'Position');
%                     delete(g);
                    g = subplot(nrows, ncols, 1);
                    
                    try % TODO there, is an error here when only some x are selected for plotting
                        hl = legend(g, Clegend, Clegendlabel, 'Location', 'SouthWest');
                        %                     lpos2 = get(hl,'Position');
                        %                     lpos2(1:2) = lpos(1:2);
                        %                     set(hl, 'Position', lpos2);
                        box(hl,'off');
                    catch
                        warning( '<arPlot2> TODO: There is an unfixed error here when only some x are selected for plotting' );
                    end
                end
                
                % optional suptitle
                if(~fastPlotTmp && exist('suptitle','file')==2 && ...
                        isfield(ar.config, 'useSuptitle') && ar.config.useSuptitle)
                    suptitle(arNameTrafo([ar.model(jm).name,': ',ar.model(jm).plot(jplot).name]))
                end
                
                % save figure
                if(saveToFile)
                    if jtype == 1
                        field = ['savePath_Fig' savepath_name{jtype} num2str(round(median(plotopt)))];
                        pfad = ['/Figures/',savepath_name{jtype},plotopt_str];
                    else
                        field = ['savePath_Fig' savepath_name{jtype}];
                        pfad = ['/Figures/',savepath_name{jtype}];
                    end
                    [ ar.model(jm).plot(jplot).(field), ...
                        ar.model(jm).plot(jplot).nRows, ...
                        ar.model(jm).plot(jplot).nCols ] = ...
                        arSaveFigure(h, ar.model(jm).plot(jplot).name, pfad);
                end
            else
                if(isfield(ar.model(jm).plot(jplot), fighandel_name{jtype}))
                    if(sum(hs==ar.model(jm).plot(jplot).(fighandel_name{jtype}))==0)
                        try %#ok<TRYNC>
                            if ( fastplot ~= 2 )
                                close(ar.model(jm).plot(jplot).(fighandel_name{jtype}))
                            end
                        end
                    end
                end
                ar.model(jm).plot(jplot).(fighandel_name{jtype}) = [];
            end
        end
        if(didPlot)
            figcount = figcount + 1;
        end
    end
end

if(nargout>0)
    varargout(1) = {hsnew};
end

