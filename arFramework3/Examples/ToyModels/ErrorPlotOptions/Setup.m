%%
arInit;
arLoadModel('model_plotErrors');
arLoadData('data_TimeCourse');
arLoadData('data_DoseResponse');
arCompileAll;

arLoadPars('IllustrationWorkspace')

%% calculate PPL
% doPPL(1,1,1,linspace(0,100,11),1);

%%
arLoad('IllustrationWorkspace');
global ar



%%
close all
dosave = true;%false;
for ploterrors = [-2,0:2]
    ar.config.ploterrors = ploterrors;
    for fiterrors = [0,1]
        ar.config.fiterrors = fiterrors;
        fprintf('fiterrors=%i, ploterrors=%i\n',ar.config.fiterrors,ar.config.ploterrors);
        arPlot2(dosave)
        drawnow
        arPlotY(dosave)
        drawnow
    end
end
