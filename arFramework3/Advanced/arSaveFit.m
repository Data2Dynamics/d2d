% function arSaveFit([name])
% 
% Saves fit history from last fit into ar.fit_hist for comparing
%
% name:     Name can be provided, otherwise it will be asked for 
%           (will be stored in ar.fit_hist(j).name)
%
% Note: ar.fit_hist will automatically appended if multiple fits will be
% saved


function arSaveFit(name)

global ar

if(~isfield(ar, 'fit_hist'))
    ar.fit_hist = struct([]);
    j = 1;
else
    j = length(ar.fit_hist) + 1;
end

if(nargin==0)
    name = ar.config.optimizers{ar.config.optimizer};
    if(ar.config.optimizer==5)
        tmpnames = arNLS;
        name = [name '_' tmpnames{ar.config.optimizerStep+1}];
    end
    addname = input(sprintf('enter name addition [%s_...]: ', name), 's');
    if(~isempty(addname))
        name = [name '_' addname];
    end
end


ar.fit_hist(j).hist = ar.fit;
ar.fit_hist(j).optimizer = ar.config.optimizer;
if(ar.config.optimizer==5)
    ar.fit_hist(j).optimizerStep = ar.config.optimizerStep;
else
    ar.fit_hist(j).optimizerStep = nan;
end
ar.fit_hist(j).config = ar.config.optim;
ar.fit_hist(j).name = name;

[~,imin] = min(ar.fit.chi2_hist + ar.fit.constr_hist);
ar.fit_hist(j).p = ar.fit.p_hist(imin,:);
