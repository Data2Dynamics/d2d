% arCheckLogXconversion(modelname,dataname,workspace,zeroVal)
% 
% 
% Example (Boehm model):
% modelname = 'FullModel';
% dataname = 'TimeCourseData';
% arCheckLogXconversion(modelname,dataname,1)
% arCheckLogXconversion(modelname,dataname,1,-10)

function arCheckLogXconversion(modelname,datanames,workspace,zeroVal)

global ar

if ~exist('modelname','var') || isempty(modelname)
    if length(ar.model)>1
        warning('Currently not supported!')
        warning('We will just return')
        return
    end
    modelname = ar.model.name;
end
if ~exist('datanames','var') || isempty(datanames)
    if length(ar.model)>1
        error('Currently not supported!')
    end
    datanames = {ar.model.data(:).name}; % does this work?
end
if ~exist('workspace','var') || isempty(workspace)
    workspace = 1;
end
if ~exist('zeroVal','var') || isempty(zeroVal)
    zeroVal = [];
end

if isstr(datanames)
    datanames = {datanames};
end

close all
clc


%% Load models & data
arInit;
arLoadModel(modelname);
for i = 1:length(datanames)
    arLoadData(datanames{i});
end
arCompileAll % ar.config.fiterrors = -1;
arFindInputs
arLoadPars(workspace)
arLin = arDeepCopy(ar);
arSave('Lin') % Only Lin workspace

% prepare logX-ing
conv = arWriteModelDef_logX(ar.model,['Models/',modelname,'_logX.def'],zeroVal);
for i = 1:length(datanames)
    arWriteDataDef_logX(ar.model.data(i),['Data/',datanames{i},'_logX.def'],conv,zeroVal);
    if exist(['Data/',datanames{i},'.xls'],'file')
        copyfile(['Data/',datanames{i},'.xls'],['Data/',datanames{i},'_logX.xls'],'f');
    elseif exist(['Data/',datanames{i},'.xlsx'],'file')
        copyfile(['Data/',datanames{i},'.xlsx'],['Data/',datanames{i},'_logX.xlsx'],'f');
    elseif exist(['Data/',datanames{i},'.csv'],'file')
        copyfile(['Data/',datanames{i},'.csv'],['Data/',datanames{i},'_logX.csv'],'f');
    else
        error('Data file %s not found.',['Data/',datanames{i}])
    end
end

% Only Log workspace
arInit
ar.config.checkForNegFluxes = false;
arLoadModel([modelname,'_logX']);
for i = 1:length(datanames)
    arLoadData([datanames{i},'_logX']);
end
arCompileAll; % ar.config.fiterrors = -1;
arFindInputs
arLoadPars(workspace)

% translate init_* to init_*LG
ind = find(~cellfun(@isempty,regexp(ar.pLabel,'init_[\w]+LG$')));
for i=1:length(ind)
    tmp = regexp(ar.pLabel(ind(i)),'init_[\w]+LG$','match');
    tref = strmatch(tmp{1}{1}(1:end-2),arLin.pLabel,'exact');
    fprintf('Translate parameter information from %s to %s...\n',arLin.pLabel{tref},ar.pLabel{ind(i)});
    ar.type(ind(i))=ar.type(tref);
    if(arLin.qLog10(tref)==1)
        ar.p(ind(i))=arLin.p(tref)/log10(exp(1));
        ar.lb(ind(i))=arLin.lb(tref)/log10(exp(1));
        ar.ub(ind(i))=arLin.ub(tref)/log10(exp(1));
        ar.mean(ind(i))=arLin.mean(tref)/log10(exp(1));
        ar.std(ind(i))=arLin.std(tref)/log10(exp(1));
        ar.qLog10(ind(i))=0;
        if(ar.type(tref)~=0)
            warning('Prior info')
        end
    else
        warning('Parameter info not convertible because of log(log)')
    end
    ar
end

arSave('Log')

% joint lin and log workspace
arInit;
ar.config.checkForNegFluxes = false;
arLoadModel(modelname);
arLoadModel([modelname,'_logX']);
for i = 1:length(datanames)
    arLoadData(datanames{i},1);
end
for i = 1:length(datanames)
    arLoadData([datanames{i},'_logX'],2);
end
arCompileAll; % ar.config.fiterrors = -1;
arFindInputs
arLoadPars(workspace)

% translate init_* to init_*LG
ind = find(~cellfun(@isempty,regexp(ar.pLabel,'init_[\w]+LG$')));
for i=1:length(ind)
    tmp = regexp(ar.pLabel(ind(i)),'init_[\w]+LG$','match');
    tref = strmatch(tmp{1}{1}(1:end-2),ar.pLabel,'exact');
    fprintf('Translate parameter information from %s to %s...\n',ar.pLabel{tref},ar.pLabel{ind(i)});
    ar.type(ind(i))=ar.type(tref);
    if(ar.qLog10(tref)==1)
        ar.p(ind(i))=ar.p(tref)/log10(exp(1));
        ar.lb(ind(i))=ar.lb(tref)/log10(exp(1));
        ar.ub(ind(i))=ar.ub(tref)/log10(exp(1));
        ar.mean(ind(i))=ar.mean(tref)/log10(exp(1));
        ar.std(ind(i))=ar.std(tref)/log10(exp(1));
        ar.qLog10(ind(i))=0;
        if(ar.type(tref)~=0)
            warning('Prior info')
        end
    else
        warning('Parameter info not convertible because of log(log)')
    end
    ar
end
arPlot
close all

arSave('BothScales')

%%
% arPrint
% arQplot('xy')
% arPlot

%% For comparison: plot some figures
arPlotLogXvsLinX


% plot(ar.res(1:32),ar.res(33:end),'.')
