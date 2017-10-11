% This function can be called after setting up the model (Setup_IL6_Full).
% Note though that it takes a long time both to compile the model as well 
% as to simulate this!

addpath PlotFunctions;
addpath Helper;

pngs = 1;
epss = 0;
tikz = 0;
addExpToFine;

if ( pngs )
    drawnow; close all;
    CoreModelPlots(1,1);
    drawnow; close all;
    APPPlot(1,1);
    drawnow; close all;
    APPTCPlots(1,1);
    drawnow; close all;
    ValidationPlot(1,1);
    drawnow; close all;
    xiaoPlot(1,1);
    drawnow; close all;
    TripleValidationPlot(1,1);
end

if ( epss )
    drawnow; close all;
    CoreModelPlots(0,1);
    drawnow; close all;
    APPPlot(0,1);
    drawnow; close all;
    APPTCPlots(0,1);
    drawnow; close all;
    ValidationPlot(0,1);
    drawnow; close all;
    xiaoPlot(0,1);
    drawnow; close all;
    TripleValidationPlot(0,1);
    drawnow; close all;
end
