% save plots

function arSavePlots(doXs, doVs, doLegends)

global ar;

if(~exist('doXs','var'))
    doXs = false;
end
if(~exist('doVs','var'))
    doVs = false;
end
if(~exist('doLegends','var'))
	doLegends = true;
end

nplots = 0;
for jm=1:length(ar.model)
    ar.model(jm).qPlotYs = false(1,length(ar.model(jm).plot));
    ar.model(jm).qPlotXs = false(1,length(ar.model(jm).plot));
    ar.model(jm).qPlotVs = false(1,length(ar.model(jm).plot));
    nplots = nplots + length(ar.model(jm).qPlotYs);
end

arWaitbar(0);
njplots = 0;
for jm=1:length(ar.model)
    for jplot=1:length(ar.model(jm).qPlotYs)
        
        njplots = njplots + 1;
        arWaitbar(njplots, nplots);
        
        ar.model(jm).qPlotYs(jplot) = 1;
        arPlot(true, false, true, false, doLegends);
        ar.model(jm).qPlotYs(jplot) = 0;
        
        if(doXs)
            ar.model(jm).qPlotXs(jplot) = 1;
            arPlot(true, false, true, false);
            ar.model(jm).qPlotXs(jplot) = 0;
        end
        
        if(doVs)
            ar.model(jm).qPlotVs(jplot) = 1;
            arPlot(true, false, true, false);
            ar.model(jm).qPlotVs(jplot) = 0;
        end
    end
end
arPlot;
arWaitbar(-1);