%% Simulate data at suggested time points
% Load Realistic model
% TimePoints have to be simulated before RealisticSimu with arRealisticDesign
% Or load a realistic model where you just want to change plog or p2 or ...

function arRealisticSimu(biomodelname,realisticname,plog,p2,dataname)

global ar
if ~exist('realisticname','var') || isempty(realisticname)
    realisticname = 'Realistic';
end
if ~exist('plog','var')
    plog = true;
end
if ~exist('p2','var')
    p2 = true;
end
if plog
    realisticname = [realisticname '_Log'];
end
if p2
    realisticname = [realisticname '_p2'];
end

%% Compile
fprintf('Compile new Observables and Parameters.\n')
arInit
ar.config.checkForNegFluxes = false;
if exist(['Models/' biomodelname],'file') || exist(['Models/',biomodelname,'.def'],'file')
    arLoadModel(biomodelname)
else
    error(['RealisticSimu.m: Couldn''t find ' biomodelname ]);
end
if exist('dataname','var') && ~isempty(dataname)
    arLoadData(dataname)
else
    arLoadData('RealisticData')
end
arCompileAll

%% Set Bounds
arSetParsScalingZero(2) %offset/scale
if plog
    arLogPars
end
if p2
    arSetParsBounds(2)
end

ar.qFit(:)=1;
ar.config.fiterrors = 1;
arLink
arSimu(false,true)

%% Simulate data on realistic t
tT = readmatrix('RealisticDesign/TimePoints.txt');
y = nan(size(tT));
for i = 1:size(tT,2) 
    arSimuData(1,1,tT(~isnan(tT(:,i)),i));
    arSimu
    y(1:size(ar.model.data.tExp,1),i) = ar.model.data.yExp(:,i);
end
[T,yExp] = artExpToVector(tT,y);
writecell([['t', ar.model.data.y];num2cell([T,yExp])],['Data' filesep 'RealisticData.xls']);

ar.model.data.tExp = T;
ar.model.data.tT = tT;
ar.model.data.yExp = yExp;
ar.model.data.yExpStd = nan(size(y));

arLink;
arSimu(true,true);
arPlot;
arSave(realisticname)
% [~,ws]=fileparts(ar.config.savepath);
% movefile(['Results' filesep ws],['Results' filesep realisticname]);
% fprintf(['Realistic simulation saved to ./Results/' realisticname '/workspace.mat \n']);