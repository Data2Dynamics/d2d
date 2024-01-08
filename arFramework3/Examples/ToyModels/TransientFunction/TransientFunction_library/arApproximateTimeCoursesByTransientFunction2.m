% fits = arApproximateTimeCoursesByTransientFunction2([doplotFolder],[boundfactor],[useRealisticTimes],[qPositive])
% 
%   RTF (retarded transient function with time-shift parameter) is fitted
%   to ODE solutions
% 
%   doplotFolder       logical indicating wether plots are made or
%                      string indicating the folder where plots are saved
%                      Default: [false]   (no plotting)
% 
%   boundfactor        a factor which can be used to increas the bounds
%                      1: like suggested for data and used in the paper 
%                      2: default
%                      4: rather broad
% 
%   useRealisticTimes   [true]
% 
%   qPositive           [false]
% 
%   This function approximates the dynamics of all model conditions and
%   states. 
% 
%   The model def is in the Examples folder of the transient function. The
%   name of the model def is TransientFunction_ForConditionFit2.def
% 
% Examples:
% fits = arApproximateTimeCoursesByTransientFunction2;
% fits = arApproximateTimeCoursesByTransientFunction2(1)
% fits = arApproximateTimeCoursesByTransientFunction2('Plots')
% 
% See also arApproximateTimeCoursesByTransientFunction

function fits = arApproximateTimeCoursesByTransientFunction2(doplot,boundfactor,useRealisticTimes,qPositive)
if ~exist('qPositive','var') || isempty(qPositive)
    qPositive = false;
end
if ~exist('useRealisticTimes','var') || isempty(useRealisticTimes)
    useRealisticTimes = true;
end

if ~exist('doplot','var') || isempty(doplot)
    doplot = false;
    plotfolder = 'plotFitTransient_with_Toffset2';
elseif ischar(doplot)
    plotfolder = doplot;
    doplot = true;
else
    plotfolder = 'plotFitTransient_with_Toffset2';
end


if ~exist('boundfactor','var') || isempty(boundfactor)
    boundfactor = 2;% In this setting (at least 20 data points, the following bounds are better)
end

plotfolder = [plotfolder,'_boundfac',num2str(boundfactor)];

global ar

if isempty(which('arFitTransient'))
    error('Please add folder to the library for the transient function, \ne.g. via %s','addpath(''...\d2d\arFramework3\Examples\ToyModels\TransientFunction\TransientFunction_library'')');
end

%% Step 1: Create data structs from conditions:
def_name = 'TransientFunction_ForConditionFit2';
def_file = [fileparts(which('arInit')),filesep,'Examples',filesep,'ToyModels',filesep,'TransientFunction',filesep,'TransientFunction_library',filesep,def_name,'.def'];
copyfile(def_file,'Models');
if ~strcmp(def_name, ar.model(end).name)
    arLoadModel('TransientFunction_ForConditionFit2');
end
fits = cell(0);
maxTF = [];
for m=1:length(ar.model)
    for ix=1:length(ar.model(m).x)
        for c=1:length(ar.model(m).condition)
            fits{end+1} = struct;
            % since each condition is fitted individually below, no condidtoin-specific prameters are required.
            fits{end}.data = arConditon2NewDataStruct2(m,c,ix,false,[],useRealisticTimes);
            fits{end}.u = ar.model(end).u;
            fits{end}.m = m;
            fits{end}.c = c;
            fits{end}.ix = ix;
            fits{end}.x = ar.model(m).x{ix};
            fits{end}.qPositive = ar.model(m).qPositiveX(ix) || sum(ar.model(m).condition(c).xFineSimu(:,ix)<0)==0;
            fits{end}.condition = struct;
            fits{end}.condition.pold = ar.model(m).condition(c).pold;
            fits{end}.condition.fp = ar.model(m).condition(c).fp;
            fits{end}.condition.p = ar.model(m).condition(c).p;
            fits{end}.label = sprintf('Model%i_%s_Condition%i',fits{end}.m,fits{end}.x,fits{end}.c);
            
            fits{end}.data.fp{4} = ['(',num2str(nanmax(fits{end}.data.tExp)),')'];  % 15.10.19, always rescale times based on tExp
        end
    end
end
fprintf('%i data-structs created. Start fitting now...\n',length(fits));
ar.model = ar.model(1:(end-1)); % remove transient function

 
%% Step2: Fit the transient functions to each data struct:
ok = cell(0);
pcatch = cell(0);

%%
% load tmp fits % if specific fits should only be fitted
%%
for d=1:length(fits)    
    try
        arInit;
        arLoadModel('TransientFunction_ForConditionFit2');
%         fits{d}.data.fp{4} = '(100)'  % for bachmann and conditions with too short tExp
        arAddDataStruct(fits{d}.data);
        arCompileAll;
        
        Initialize_FitTransient2(boundfactor,[],qPositive);
        ar.fit_transient.doReduction = false;
        arFitTransient

        % 2 loops of LHS: First 10, then 20 fits are done
        SDest = ar.p(ar.fit_transient.indp.sd);
        if ar.qLog10(ar.fit_transient.indp.sd)==1
            SDest = 10.^SDest;
        end
        if sum(~isnan(ar.model.data.yExp))>1 && range(ar.model.data.yExp)/10<SDest  && range(ar.model.data.yExp)>1e-10
            fprintf('Fit %i out of %i has large SD: LHS(10) started ...\n',d,length(fits));
            arFitLHS(10)  % 10 fits should be enough for 7 parameters!
            
            arFits(ar.ps); % I don't know why this works (so well)
                
            % if not yet good, make arFitLHS(20)
            SDest = ar.p(ar.fit_transient.indp.sd);
            if ar.qLog10(ar.fit_transient.indp.sd)==1
                SDest = 10.^SDest;
            end            
            if sum(~isnan(ar.model.data.yExp))>1 && range(ar.model.data.yExp)/10<SDest  && range(ar.model.data.yExp)>1e-10
                fprintf('Fit %i out of %i has still large SD: LHS(20) started ...\n',d,length(fits));
                arFitLHS(20)
            end
        end
        
        arSimu(false,true,true);
        arCalcMerit
        
        fits{d}.p = ar.p;
        indLog = find(ar.qLog10==1);
        fits{d}.p(indLog) = 10.^fits{d}.p(indLog); 
        
        fits{d}.merit = arGetMerit;
        fits{d}.pLabel = ar.pLabel;
        fits{d}.tFine = ar.model.data.tFine;
        fits{d}.yFineSimu = ar.model.data.yFineSimu;
        fits{d}.ystdFineSimu = ar.model.data.ystdFineSimu;
        
        fprintf('Fit %4i out of %5i done.\n',d,length(fits));
        if rem(d,10)==0 % save every tenth times
            save fits_tmp fits
        end
        ok{end+1} = struct('p',ar.p,'lb',ar.lb,'ub',ar.ub,'yExp',ar.model.data.yExp,'tExp',ar.model.data.tExp);
        
        if doplot        
            plotFitTransient(fits(d),plotfolder);
        end

    catch ERR
        pcatch{end+1} = struct('p',ar.p,'lb',ar.lb,'ub',ar.ub);
        save error ok pcatch
    end
end
save arApproximateTimeCoursesByTransientFunction2
arPlot


%% calculate the approximation error from the fitted SD
for i=1:length(fits)
    if(isfield(fits{i},''))
        fits{i}.approxErr = 10.^fits{i}.p(4)./range(fits{i}.data.yExp);
    else
        fits{i}.approxErr = NaN*range(fits{i}.data.yExp);
    end
end


