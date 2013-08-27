function arShowAllFits

global ar

close all

for jm=1:length(ar.model)
    ar.model(jm).qPlotYs(:) = 0;
    ar.model(jm).qPlotXs(:) = 0;
    ar.model(jm).qPlotVs(:) = 0;
end
arPlotY;

for jm=1:length(ar.model)
    for jplot=1:length(ar.model(jm).plot)
        ar.model(jm).qPlotYs(jplot) = 1;
        arPlotY;
        waitforbuttonpress;
        ar.model(jm).qPlotYs(jplot) = 0;
    end
end

close all