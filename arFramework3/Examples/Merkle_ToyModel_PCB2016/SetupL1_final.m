arInit;

arLoadModel('model_template');

arLoadData('data_template1');
arLoadData('data_template2');

arCompileAll;

ar.config.useNewPlots=1;

ar.p(1) = 0;
ar.p(2) = log10(.2);
ar.p(3) = log10(.5);
ar.p(4) = log10(.5);
ar.p(5) = log10(2);
ar.p(6:7) = 0;
ar.p(8) = log10(.2);

l1pars = 5:8;

p_init = ar.p;

ar.model.data(2).yExp(:) = nan;
ar.model.data(2).ndata(:) = 0;
ar.model.data(1).tLim = ar.model.tLim;
ar.model.data(2).tLim = ar.model.tLim;
arLink

ar.model.data(1).y = {'Protein','ppProtein'};
ar.model.data(2).y = {'Protein','ppProtein'};
ar.model.data(1).yNames = {};
ar.model.data(2).yNames = {};
ar.model.x = {'Protein','pProtein','ppProtein'};

rng(0); % Make results reproducible
ar.p = p_init;
ar.qFit(:) = 1;

arSimuData([],[],linspace(0,60,21))

arFit
arPlot

removeInds = 1:2:21;

ar.model.data(1).yExp(removeInds,2) = NaN;
ar.model.data(2).yExp(removeInds,2) = NaN;

arFit
arPlot
arPrint

p_fit = ar.p;
ar.p = p_init;
ar.model.qPlotYs(:) = 1;
ar.model.qPlotXs(:) = 1;
arPlot
close all
ar.model.qPlotXs(:) = 0;
arPlot

ar.p = p_fit;


l1Init(l1pars)

linv = logspace(-4,2,20);
linv = linv(end:-1:1);

l1Scan([],linv)

l1Unpen
l1SelectOpt
l1Plot

arPlot