% arFitLHS([n], [randomseed], [log_fit_history], [backup_save], [use_cluster])
%
% Multi-Start optimization fit sequence using random sampling from prios 
% or latin hyper cube sampling 
%
%   n                       number of fits	[10]
%   randomseed              random seed used for initial guesses
%                           which enters rng(randomseed), if empty a random
%                           random seed it used (i.e. rng('shuffling'))
%                           []  (if left empty, seed is internally chosen)
%   log_fit_history         if true, details of all individuals fits from ar.fit
%                           are stored in ar.fit_hist 
%                           [false]
%   backup_save             if true, ar struct with the fitting sequence is 
%                           stored in arFits_backup.mat
%                           [false]
%   use_cluster             option to use the function on Matlab Distributed 
%                           Compute Server (MDCS), but *not* on bwHPC Custers
%                           [false]
%
%   ar.config.useLHS        if true, latin hyper cube (LHS) sampling is used 
%                           for drawing the inital guesses. Otherwise,
%                           random samples are drawn randomly from the
%                           priors or box contraints
%                           [false]
%   ar.config.restartLHS    if set to 1 and integration is not feasible, the 
%                           fit is restarted with a new random initial guess.
%                           If set to 0, non-feasible fits are possible.
%                           [0]
%
% This function draws initial guesses from the parameter space, either 
% using 'latin hyper cube sampling' (LHS) or 'random sampling' for a given 
% number of multi start fits. A random seed can either be choosen directely,
% or it is set internally if argument is left empty. Afterwards, the fit 
% sequence is started using arFit and the best fit of this sequence is used 
% as new best fit (i.e. in ar.p). Results of the seuqence are available in 
% ar.ps, ar.ps_start, ar.chi2s, ar.chi2s_start, etc.. If log_fit_history is 
% set to 'true', detailed results of each fit are logged and stored in 
% ar.fit_hist  Results (i.e. waterfall plots) can be plotted afterwards 
% using arPlotFits or arPlotChi2s. Additionally results from 
% arLhsSampleSizeCalculation are shown after the fit sequence if more than 
% 50 fits were initiated.
%
% Examples:
%
%     arFitLHS(25)
% 
%     25 fits are performed (initil guesses drawn randomly with internally set 
%     random seed).
% 
%     ar.config.useLHS = 0;
%     arFitLHS(100,1337,[],1)
%     arPlotChi2s;
% 
%     100 fits from latin hyper cube sampling are drawn with the random seed 
%     '1337' and the resulting ar struct is stored in arFits_backup.mat. 
%     Afterwards, the waterfall plot is 
% 
% See also arFit, arLhsSampleSizeCalculation, arPlotFits, arPlotChi2s, arChi2LHS


function arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)

if(~exist('n','var'))
    n = 10;
end
if(~exist('randomseed','var') || isempty(randomseed))
    randomseed = [];
end
if(~exist('log_fit_history','var') || isempty(log_fit_history))
    log_fit_history = false;
end
if(~exist('backup_save','var') || isempty(backup_save))
    backup_save = false;
end
if(~exist('use_cluster','var') || isempty(use_cluster))
    use_cluster = false;
end

global ar
if ~isfield(ar.config,'restartLHS')
    ar.config.restartLHS = 0;
end
if(~isfield(ar.config,'nCVRestart') || isnan(ar.config.nCVRestart))
    nCV_bkp = NaN;
    ar.config.nCVRestart = 1;
else
    nCV_bkp = ar.config.nCVRestart;
end

% generate random values
ps = arRandomPars(n, randomseed);

if(~use_cluster)
    arFits(ps, log_fit_history, backup_save);
else
    arFitsCluster(ps, log_fit_history, backup_save);
end
ar.config.nCVRestart = nCV_bkp;
if ar.config.restartLHS ==1 && isempty(randomseed)
    indnan = find(isnan(ar.chi2s));
    counter = 1;
    ar.lhsRepeats = ones(size(ar.chi2s));
    while ~isempty(indnan) && counter<10
        counter = counter+1;
    
        fprintf('Repeat fits %i\n', indnan);

        ar.lhsRepeats(indnan) = ar.lhsRepeats(indnan) +1;
        
        pstmp = NaN(size(ps));
        psNeu = arRandomPars(n, randomseed);
        pstmp(indnan,:) = psNeu(indnan,:);
        
        if(~use_cluster)
            arFits(pstmp, log_fit_history, backup_save);
        else
            arFitsCluster(pstmp, log_fit_history, backup_save);
        end
        indnan = find(isnan(ar.chi2s));
    end
else
    if(sum(isnan(ar.chi2s))>0)
        fprintf('\nSome fits yield parameters where ODE intergration was not feasible.\n');
        fprintf('Such fits can be automatically restarted with another intial guess by setting\n');
        fprintf('   ar.config.restartLHS = 1\n\n')
    end
end

try
    ar.LhsSampleSizeCalculation = arLhsSampleSizeCalculation;
catch
    disp(lasterror)
end

