% hs = arPlot([saveToFile], [fastPlot], [silent], [sensi], [doLegends], [dynamics], [hs])
%
% Plot models and datasets
%
% saveToFile    Save plot to file               [false]
% fastPlot      Do fast plotting                [false]
% silent        Plot without printing text      [false]
% sensi         Uses sensitivities              [false] (replace position of deprecated evalfun)
% doLegends     Print the legends               [true]
% dynamics      Simulate dynamics               [true]
% hs            Plot to custom figure handles   []
%
% ar.config.fiterror     0: Data is plotted as fitted (default).
%                       -1: Plot prediction bands as calculated by PPL.
%                       -2: Do not plot errors
% ar.config.fiterrors   -1: Only experimental error bars used for fitting.
%                        0: Use experimental errors by default and revert
%                           to the error model for data points that have no
%                           experimental error specified (default).
%                        1: Only the error model is used for fitting.
%                           Experimental errors specified in the data sheet
%                           are ignored.
% ar.config.ploterrors   0: Observables are plotted as fitted (default).
%                        1: Data uncertainty is plotted as error bar.
%                        2: Only error bands are plotted.
%
% By default arPlot simulates the model without sensitivities and subsequently 
% plots all the enabled observables, states/derived variables and
% fluxes. Which conditions are plotted can be set with arPlotter. 
%
% Note: arPlot2 also plots model simulations using different rendering code.
% Note: For observables, axis transformations can be set in 
%       ar.model(jm).plot(jplot).xtrafo and ar.model(jm).plot(jplot).ytrafo
%       by specifying single input/output anonymous functions here.
%
% See also arPlotter, arPlot2, arPlotY, arPlotX, arPlotV.

function varargout = arPlot(saveToFile, fastPlot, silent, sensi, doLegends, dynamics, hs)

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
        arSimu(sensi, true, dynamics);
    catch err_id
        if(~silent)
            disp(err_id.message);
        end
    end
    try 
        if(silent)
            arCalcMerit(sensi, ar.p(ar.qFit==1), dynamics);
        else
            arCalcMerit(sensi);
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