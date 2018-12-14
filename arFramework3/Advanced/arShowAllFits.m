% arShowAllFits()
%
% Show successively all plots. Next plot occurs after a key press.
% After the function call all qPlotYs, qPlotXs and qPlotVs are zero

function arShowAllFits

global ar

close all

for jm=1:length(ar.model)
    ar.model(jm).qPlotYs(:) = 0;
    ar.model(jm).qPlotXs(:) = 0;
    ar.model(jm).qPlotVs(:) = 0;
end

for jm=1:length(ar.model)
    for jplot=1:length(ar.model(jm).plot)
        ar.model(jm).qPlotYs(jplot) = 1;
        if(isfield(ar.config, 'useNewPlots') && ar.config.useNewPlots)
            arPlot2;
        else
            arPlotY;
        end
        waitforbuttonpress;
        ar.model(jm).qPlotYs(jplot) = 0;
    end
end

close all