% hs = arPlot([saveToFile], [fastPlot], [silent], [evalfun], [doLegends], [dynamics], [hs])
%
% Plot models and datasets
%
% saveToFile    Save plot to file               [false]
% fastPlot      Do fast plotting                [false]
% silent        Plot without printing text      [false]
% evalfun       Evaluate function               [true]  (DEPRECATED)
% doLegends     Print the legends               [true]
% dynamics      Simulate dynamics               [true]
% hs            Plot to custom figure handles   []
%
% arPlot simulates the model without sensitivities and subsequently 
% plots all the enabled observables, states/derived variables and
% fluxes. Which observables/states and fluxes are plotted can be set 
% with arPlotter. 
%
% Note: arPlot2 also plots model simulations using different rendering code.
% Note: For observables, axis transformations can be set in 
%       ar.model(jm).plot(jplot).xtrafo and ar.model(jm).plot(jplot).ytrafo
%       by specifying single input/output anonymous functions here.
%
% See also arPlotter, arPlot2, arPlotY, arPlotX, arPlotV.

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
evalfun = true;
% Always potentially evaluate the model. This will not hurt performance since 
% arSimu checks against the cache. Running with false will however cause problems.
% if(~exist('evalfun','var'))
%     evalfun = true;
% end
if(~exist('doLegends','var'))
    if(isfield(ar.config, 'showLegends'))
        doLegends = ar.config.showLegends;
    else
        doLegends = true;
    end
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
            arCalcMerit(false, ar.p(ar.qFit==1), dynamics);
        else
            arCalcMerit(false);
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