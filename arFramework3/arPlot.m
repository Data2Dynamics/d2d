% Plot models and datasets
%
% hs = arPlot(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics, hs)
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

function varargout = arPlot(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics, hs)

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


if(isfield(ar.config, 'useNewPlots') && ar.config.useNewPlots)
    hs = arPlot2(saveToFile, fastPlot, silent, evalfun, doLegends, dynamics, hs);
    if(nargout>0)
        varargout(1) = {hs};
    end
    return;
else
    if(nargout>0)
        varargout(1) = {hs};
    end
end

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
            arChi2(false, ar.p(ar.qFit==1), dynamics);
        else
            arChi2;
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

arPlotY(saveToFile, fastPlot, doLegends);
arPlotX(saveToFile, fastPlot);
arPlotV(saveToFile, fastPlot);