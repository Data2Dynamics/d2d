% fits = arFitTransientFunction2(dat,[doplotFolder],[boundfactor])
%
%   Transient function with time-shift parameter.
%   This function fits experimental data provided as arguments dat.
%
%   dat.tExp
%   dat.yExp
%   dat.ystd
%
%   doplotFolder       logical indicating wether plots are made or
%                      string indicating the folder where plots are saved
%                      Default: [false]   (no plotting)
%
%   boundfactor        [1]
%                      Factor which might be used to enlarge the default
%                      bounds
%
% Example:
%
% dat.tExp = ar.model(1).data(1).tExp;
% dat.yExp = ar.model(1).data(1).tExp;
% dat.ystd = ar.model(1).data(1).tExpStd;
% fits = arFitTransientFunction2(dat)

function res = arFitTransientFunction2(dat,doplot,boundfactor,qPositive)
if ~exist('qPositive','var') || isempty(qPositive)
    qPositive = false;  % Is the truth known to be positive?
end
if ~exist('doplot','var') || isempty(doplot)
    doplot = false;
    plotfolder = 'plotFitTransient_with_Toffset';
elseif ischar(doplot)
    plotfolder = doplot;
    doplot = true;
else
    plotfolder = 'plotFitTransient_with_Toffset';
end

if ~exist('boundfactor','var') || isempty(boundfactor)
    boundfactor = 1;% In this setting (at least 20 data points, the following bounds are better)
end

% plotfolder = [plotfolder,'_boundfac',num2str(boundfactor)];

if isempty(which('arFitTransient'))
    error('Please add folder to the library for the transient function, \ne.g. via %s','addpath(''...\d2d\arFramework3\Examples\ToyModels\TransientFunction\TransientFunction_library'')');
end


%% Step 1: Create data structs from conditions:
%def_file = [fileparts(which('arInit')),filesep,'Examples',filesep,'ToyModels',filesep,'TransientFunction',filesep,'TransientFunction_library',filesep,'TransientFunction_ForConditionFit2.def'];
def_file = [fileparts(mfilename('fullpath')),filesep,'TransientFunction.def'];
copyfile(def_file,'Models','f');

%% Check for inf
if any(isinf(dat.yExp))
    dat.yExp(isinf(dat.yExp) & dat.yExp>0) = max(dat.yExp(~isinf(dat.yExp)))*2;
    dat.yExp(isinf(dat.yExp) & dat.yExp<0) = min(dat.yExp(~isinf(dat.yExp)))*2;
    warning('Set inf value to 2* max/min (depending on sign).')
end
if any(isinf(dat.yExp))
    dat.ystd(isinf(dat.ystd)) = nan;
    warning('Set inf measurement value to nan.')
end

global ar
arInit
arLoadModel('TransientFunction');

args = cell(0);
args{end+1} = 'tExp';       args{end+1} = dat.tExp;
args{end+1} = 'yExp';       args{end+1} = dat.yExp;
args{end+1} = 'yExpStd';       args{end+1} = dat.ystd; %*NaN;  % sd not used
args{end+1} = 'yNames';          args{end+1} = ar.model.yNames;
args{end+1} = 'y';          args{end+1} = ar.model.y;

maxt_TF = range(dat.tExp);
D = arCreateDataStruct(1,{'maxt_TF'},{['(',num2str(maxt_TF),')']},args{:});
arAddDataStruct(D);
arCompileAll;


Initialize_FitTransient2(boundfactor,[],qPositive);
%         if D.qPositive==1
% if qPositive==1
%     ar.lb(ar.fit_transient.indp.offset) = 0;
%     ar.p(ar.fit_transient.indp.offset) = max(ar.p(ar.fit_transient.indp.offset),ar.lb(ar.fit_transient.indp.offset));
%     ar.ub(ar.fit_transient.indp.offset) = max(0,ar.ub(ar.fit_transient.indp.offset));
%     
%     ar.fit_transient.bounds.lb(ar.fit_transient.indp.offset) = max(0,ar.fit_transient.bounds.lb(ar.fit_transient.indp.offset));
%     ar.fit_transient.boundsNeg.lb(ar.fit_transient.indp.offset) = max(0,ar.fit_transient.boundsNeg.lb(ar.fit_transient.indp.offset));
% end
arFitTransient
SDest = ar.p(ar.fit_transient.indp.sd);
if sum(ar.qFit==1)> 1 && ar.qLog10(ar.fit_transient.indp.sd)==1
    SDest = 10.^SDest;
end
if sum(ar.qFit==1)>1 && sum(~isnan(ar.model.data.yExp))>1 && range(ar.model.data.yExp)/10<SDest  && range(ar.model.data.yExp)>1e-10
    arFitLHS(10)  % 10 fits should be enough for 7 parameters!
end

arSimu(false,true,true);
arCalcMerit

res = struct;
res.data = dat;
res.p = ar.p;
res.merit = arGetMerit;
res.pLabel = ar.pLabel;
res.tFine = ar.model.data.tFine;
res.yFineSimu = ar.model.data.yFineSimu;
res.ystdFineSimu = ar.model.data.ystdFineSimu;
try
    res.label = dat.label;
catch
    res.label = '';
end

res.pRescaled = ar.p; % time parameters rescaled to original time scale of the data
ind = find(contains(ar.pLabel,{'timescale_','toffset_'}));
for i=1:length(ind)
    if ar.qLog10(ind(i))
        res.pRescaled(ind(i)) = log10(maxt_TF* 10^res.pRescaled(ind(i))/10 );
    else
        res.pRescaled(ind(i)) = maxt_TF* res.pRescaled(ind(i))/10 ;    
    end
end

%         ok{end+1} = struct('p',ar.p,'lb',ar.lb,'ub',ar.ub,'yExp',ar.model.data.yExp,'tExp',ar.model.data.tExp);

if doplot
    plotFitTransient(res,plotfolder);
end

% calculate approximation error
res.approxErr = 10.^res.p(4)./range(res.data.yExp);


