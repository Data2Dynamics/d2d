% arCheckLogXconversion(modelname,dataname,workspace,zeroVal)
% 
% 
% Example (Boehm model):
% modelname = 'FullModel';
% dataname = 'TimeCourseData';
% arCheckLogXconversion(modelname,dataname,1)
% arCheckLogXconversion(modelname,dataname,1,-10)

function arCheckLogXconversion(modelname,dataname,workspace,zeroVal)
if ~exist('zeroVal','var') || isempty(zeroVal)
    zeroVal = [];
end
if ~exist('workspace','var') || isempty(workspace)
    workspace = 1;
end

global ar

close all
clc

arInit;
arLoadModel(modelname);
arLoadData(dataname,1);
arCompileAll

conv = arWriteModelDef_logX(ar.model,['Models/',modelname,'_logX.def'],zeroVal);
arWriteDataDef_logX(ar.model.data,['Data/',dataname,'_logX.def'],conv,zeroVal);

if exist(['Data/',dataname,'.xls'],'file')
    copyfile(['Data/',dataname,'.xls'],['Data/',dataname,'_logX.xls'],'f');
elseif exist(['Data/',dataname,'.csv'],'file')
    copyfile(['Data/',dataname,'.csv'],['Data/',dataname,'_logX.csv'],'f');
else
    error('Data file %s not found.',['Data/',dataname])
end

% Load models & data
arInit;
ar.config.checkForNegFluxes = false;
% arLoadModel('ModelImmuneCells');
arLoadModel(modelname);
arLoadModel([modelname,'_logX']);
arLoadData(dataname,1);
arLoadData([dataname,'_logX'],2);
arCompileAll;

% ar.config.fiterrors = -1;

arLoadPars(workspace)

%% translate init_* to init_*LG
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

close all
%%
arPrint
arQplot('xy')
arPlot

%% Comparison
[~,ia,ib] = intersect(ar.model(1).x,ar.model(2).z,'stable');
maxdiff = [];
maxreldiff = [];
for i=1:length(ia)
    tFine = ar.model(1).condition(1).tFine;
    xlog = ar.model(1).condition(1).xFineSimu(:,ia(i));
    zlog = interp1(ar.model(2).condition(1).tFine,ar.model(2).condition(1).zFineSimu(:,ib(i)),tFine);
    loglog(xlog,...
        zlog,...
        '.')
    hold on
    maxdiff(i) = max(abs(xlog-zlog));
    maxreldiff(i) = max(abs(xlog-zlog)./((xlog+zlog)/2));
    set(gca,'ColorOrderIndex',get(gca,'ColorOrderIndex')+1);
end
maxdiff
maxreldiff
abplot(1,0);
axis tight
xlabel('standard x')
ylabel('log(x)');
legend(ar.model(1).x{:},'Location','SouthEast');




% plot(ar.res(1:32),ar.res(33:end),'.')
