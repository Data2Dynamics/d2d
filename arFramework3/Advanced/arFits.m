% arFits(ps, [log_fit_history], [backup_save], [prefunc], [postfunc])
%
% Performing a multistart fit sequence with given initial guesses using arFit
%
%   ps                  matrix of parameter values for fitting sequence
%                       [nfit x np]
%   log_fit_history     if true, details of all individuals fits from 
%                       ar.fit are stored in ar.fit_hist 
%                       [false]
%   backup_save         if true, ar struct with the fitting sequence is 
%                       stored in arFits_backup.mat
%                       [false]
%   prefunc             function to be called before each fitting
%   postfunc            function to be called after each fitting
% 
% Performs a fitting sequence of the provided matrix of initial values for 
% the parameter values using arFit. ps has to have the right dimensions 
% (rows: number of fits, columns: number of fitted parameters).
% if 
%   1) ps has the same size as ar.ps or is larger  AND
%   2) contains rows with only NaN, then this fits corresponding to NaN
%   rows are not performed and the old fit result is maintained. 
%   This enables overwriting fits where integration was not feasible (e.g. in arFitLHS)
%   and adding new fits.
% 
%
% Examples:
% % 1) Standard call:
% ps = ones(100,1)*ar.p;
% ps(:,1) = linspace(-1,1,100);
% arFits(ps)
% 
% % 2) Add new Fits to result of 1)
% ps2 = ps;
% ps2(:,2) = linspace(-1,1,100);
% arFits([NaN(size(ar.ps));ps2])
% 
% % 3) Restart a single Fit (e.g. the 5th) with altered tolerances:
% ar.config.atol = ar.config.atol/10;
% ps3 = NaN(size(ar.ps));
% ps3(5,:) = ar.ps_start(5,:);
% arFits(ps3)
% 
% See also arFit, arFitLHS, arChi2LHS

function arFits(ps, log_fit_history, backup_save, prefunc, postfunc)

global ar

if(~exist('log_fit_history','var') || isempty(log_fit_history))
    log_fit_history = false;
end
if(~exist('backup_save','var') || isempty(backup_save))
    backup_save = false;
end
if(~exist('prefunc','var') || isempty(prefunc))
    prefunc = [];
end
if(~exist('postfunc','var') || isempty(postfunc))
    postfunc = [];
end

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

dop = find(sum(~isnan(ps),2)==size(ps,2));

n = length(dop);
if ~isfield(ar,'ps') || size(ps,1)<size(ar.ps,1)
    replaceOld = true;
else
    replaceOld = false;
end

if replaceOld
    ar.ps_start = ps;
    ar.ps = nan(size(ps));
    ar.ps_errors = nan(size(ps));
    ar.chi2s_start = nan(1,size(ps,1));
    ar.chi2sconstr_start = nan(1,size(ps,1));
    ar.chi2s = nan(1,size(ps,1));
    ar.chi2sconstr = nan(1,size(ps,1));
    ar.exitflag = nan(1,size(ps,1));
    ar.timing = nan(1,size(ps,1));
    ar.fun_evals = nan(1,size(ps,1));
    ar.iter = nan(1,size(ps,1));
    ar.optim_crit = nan(1,size(ps,1));
    if(isfield(ar.config,'logFitting') && ar.config.logFitting) 
        if(ar.config.logFitting)
            ar.optimLogs = cell(1,size(ps,1));
        end
    end
else
    ar.ps_start(dop,:) = ps(dop,:);
    ar.ps(dop,:) = nan(size(ps(dop,:)));
    ar.ps_errors(dop,:) = nan(size(ps(dop,:)));
    ar.chi2s_start(dop) = nan(1,size(ps(dop,:),1));
    ar.chi2sconstr_start(dop) = nan(1,size(ps(dop,:),1));
    ar.chi2s(dop) = nan(1,size(ps(dop,:),1));
    ar.chi2sconstr(dop) = nan(1,size(ps(dop,:),1));
    ar.exitflag(dop) = nan(1,size(ps(dop,:),1));
    ar.timing(dop) = nan(1,size(ps(dop,:),1));
    ar.fun_evals(dop) = nan(1,size(ps(dop,:),1));
    ar.iter(dop) = nan(1,size(ps(dop,:),1));
    ar.optim_crit(dop) = nan(1,size(ps(dop,:),1));
    if(isfield(ar.config,'logFitting') && ar.config.logFitting) 
        if(~isfield(ar,'optimLogs'))
            ar.optimLogs = cell(1,n);
        else
            ar.optimLogs(dop) = cell(1,length(dop));
        end
    end

end

pReset = ar.p;
try
    arCalcMerit(true,ar.p(ar.qFit==1));
    chi2Reset = arGetMerit('chi2fit') + arGetMerit('chi2constr');
catch ERR
    if strcmp(ERR.identifier,'d2d:arCollectRes:NaN_in_res')~=1 && strcmp(ERR.identifier,'d2d:arCollectRes:NaN_in_sres')~=1
        rethrow(ERR)
    else % do not throw an error if NaNs are in residuals for the initial guess
        chi2Reset = Inf;
        warning(ERR.message)
    end
end

if(log_fit_history)
    ar.fit_hist = [];
end


arWaitbar(0);

global handles
handles.stop_now = 0;
buttonStop

for j=1:n
    arWaitbar(j, n);
    if j > 1
        if isfield(ar.config,'showLiveWaterfall') && ar.config.showLiveWaterfall ==1
            arPlotChi2s
            fig = gcf;
            fig.Children(3).XLim = [0 length(ar.ps)];
        end 
    end
    
    if handles.stop_now==1  
        arFprintf(1,'The fitting was stopped after %i/%i iterations manually by the user.',j-1,n);
        break;
    end
    
    ar.p = ps(dop(j),:);
    if(isfield(ar.config,'useDouble') && ar.config.useDouble==1)
        ar.p(ar.iref) = ar.p(ar.iprimary);
    end
    
    tic;
    try
        arCalcMerit(true,ar.p(ar.qFit==1));
        ar.chi2s_start(dop(j)) = arGetMerit('chi2fit');
        ar.chi2sconstr_start(dop(j)) = arGetMerit('chi2constr');
        if ~isempty(prefunc)
            try
                feval( prefunc );
            catch
                arFprintf(1, 'Error: Failure calling pre-fitting function');
            end
        end
        arFit(true);
        if ~isempty(postfunc)
            try
                feval( postfunc );
            catch
                arFprintf(1, 'Error: Failure calling post-fitting function');
            end
        end       
        ar.ps(dop(j),:) = ar.p;
        ar.chi2s(dop(j)) = arGetMerit('chi2fit');
        ar.chi2sconstr(dop(j)) = arGetMerit('chi2constr');
        ar.exitflag(dop(j)) = ar.fit.exitflag;
        ar.fun_evals(dop(j)) = ar.fit.fevals;
        ar.iter(dop(j)) = ar.fit.iter;
        ar.optim_crit(dop(j)) = ar.firstorderopt;
    catch exception
        ar.chi2s(dop(j)) = inf;
        ar.ps_errors(dop(j),:) = ar.p;
        fprintf('fit #%i: %s\n', dop(j), exception.message);
    end

    ar.timing(dop(j)) = toc;
    if(isfield(ar, 'fit'))
        if(isfield(ar.fit,'optimLog'))  % coincides with ar.config.logFitting
            ar.optimLogs{dop(j)} = ar.fit.optimLog;
        end
    end
        
    if(log_fit_history)
        name = ar.config.optimizers{ar.config.optimizer};
        if(ar.config.optimizer==5)
            tmpnames = arNLS;
            name = [name '_' tmpnames{ar.config.optimizerStep+1}]; %#ok<AGROW>
        end
        
        ar.fit_hist(dop(j)).hist = ar.fit;
        ar.fit_hist(dop(j)).optimizer = ar.config.optimizer;
        if(ar.config.optimizer==5)
            ar.fit_hist(dop(j)).optimizerStep = ar.config.optimizerStep;
        else
            ar.fit_hist(dop(j)).optimizerStep = nan;
        end
        ar.fit_hist(dop(j)).config = ar.config.optim;
        ar.fit_hist(dop(j)).name = [name '_' sprintf('run%i', dop(j))];
        
        [~,imin] = min(ar.fit.chi2_hist + ar.fit.constr_hist);
        ar.fit_hist(dop(j)).p = ar.fit.p_hist(imin,:);
    end
    if(backup_save)
        save('arFits_backup.mat','ar');
    end    
end

arFprintf([],'total fitting time: %s\n', secToHMS(sum(ar.timing(~isnan(ar.timing)))));
arFprintf([],'mean fitting time: %s\n', secToHMS(10^mean(log10(ar.timing(~isnan(ar.timing))))));
arWaitbar(-1);

if(chi2Reset>min(ar.chi2s + ar.chi2sconstr))
    [chi2min,imin] = min(ar.chi2s + ar.chi2sconstr);
    ar.p = ar.ps(imin,:);
    if ar.config.fiterrors == -1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)==0) % if no error parameters fitted
        arFprintf([],'selected best fit #%i with %f (old = %f)\n', imin, 2*(ar.ndata+ar.nconstr)*log(sqrt(2*pi)) + chi2min, 2*(ar.ndata+ar.nconstr)*log(sqrt(2*pi)) + chi2Reset);
    else
        arFprintf([],'selected best fit #%i with %f (old = %f)\n', ...
            imin, 2*ar.ndata*log(sqrt(2*pi)) + chi2min, 2*ar.ndata*log(sqrt(2*pi)) + chi2Reset);
    end
else
    arFprintf([],'did not find better fit\n');
    ar.p = pReset;
end
try
    arCalcMerit(true,[]);
catch ERR
    if strcmp(ERR.identifier,'d2d:arCollectRes:NaN_in_res')~=1 && strcmp(ERR.identifier,'d2d:arCollectRes:NaN_in_sres')~=1
        rethrow(ERR)
    else % do not throw an error if NaNs are in residuals for the initial guess
        warning(ERR.message)
    end
end


function buttonStop

persistent batchmode
if isempty(batchmode)
    ss = get(0, 'ScreenSize');
    if max(ss)<2
        batchmode = true;
    else
        batchmode = false;
    end
end

if batchmode 
    return % do not show Waitbar in batch mode
end


% Create a figure window
fig = uifigure;
fig.Position = [608   258   321   118];

% Create a UI axes
% ax = uiaxes('Parent',fig,...
%             'Units','pixels',...
%             'Position', [104, 123, 300, 201]);   

% Create a push button
btn = uibutton(fig,'push',...
               'Position',[100 50 100 22],...
               'ButtonPushedFcn', @(btn,event) stopButtonPushed(fig));
btn.Text = 'Stop fits';


% Create the function for the ButtonPushedFcn callback
function stopButtonPushed(fig2)
       global handles
       handles.stop_now = 1;
       close(fig2)
