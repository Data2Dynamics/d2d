% fit sequence using 
%   - latin hyper cube sampling (ar.config.useLHS = true)
%   - random sampling from prior
%
% arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)
%
% n:                number of runs      [10]
% randomseed:                           rng(randomseed)
% log_fit_history                       [false]
% backup_save                           [false]
% use_cluster                           [false]
% 
%   ar.config.restartLHS = 0: Default, non-feasible fits are possible.
%   ar.config.restartLHS = 1: If integration is not feasible, the fit is
%                             restarted with a new random initial guess.

function arFitLHS(n, randomseed, log_fit_history, backup_save, use_cluster)

if(~exist('n','var'))
    n = 10;
end
if(~exist('randomseed','var'))
    randomseed = [];
end
if(~exist('log_fit_history','var'))
    log_fit_history = false;
end
if(~exist('backup_save','var'))
    backup_save = false;
end
if(~exist('use_cluster','var'))
    use_cluster = false;
end

global ar
if ~isfield(ar.config,'restartLHS')
    ar.config.restartLHS = 0;
end

% generate random values
ps = arRandomPars(n, randomseed);

if(~use_cluster)
    arFits(ps, log_fit_history, backup_save);
else
    arFitsCluster(ps, log_fit_history);
end

if ar.config.restartLHS ==1 && isempty(randomseed)
    indnan = find(isnan(ar.chi2s));
    counter = 1;
    ar.lhsRepeats = ones(size(ar.chi2s));
    while ~isempty(indnan) && counter<10
        counter = counter+1;
    
        fprintf('Repeat fits');
        fprintf('%i, ',indnan);
        fprintf('\n');

        ar.lhsRepeats(indnan) = ar.lhsRepeats(indnan) +1;
        
        pstmp = NaN(size(ps));
        psNeu = arRandomPars(n, randomseed);
        pstmp(indnan,:) = psNeu(indnan,:);
        
        if(~use_cluster)
            arFits(pstmp, log_fit_history, backup_save);
        else
            arFitsCluster(pstmp, log_fit_history);
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


