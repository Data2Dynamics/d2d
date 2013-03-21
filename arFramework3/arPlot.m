% Plot models and datasets
%
% arPlot(saveToFile, fastPlot, silent, evalfun, doLegends)
%
% saveToFile    [false]
% fastPlot      [false]
% silent        [false]
% evalfun       [false]
% doLegends      [true]

function arPlot(saveToFile, fastPlot, silent, evalfun, doLegends)

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

if(evalfun)
    try
        if(silent)
            arChi2(false);
        else
            arChi2;
        end
    catch error_id
        disp(error_id.message);
    end
    try
        arSimu(false, true);
    catch error_id
        disp(error_id.message);
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

