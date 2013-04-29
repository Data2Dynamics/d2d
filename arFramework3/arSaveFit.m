function arSaveFit

global ar

if(~isfield(ar, 'fit_hist'))
    ar.fit_hist = struct([]);
    j = 1;
else
    j = length(ar.fit_hist) + 1;
end

name = input(sprintf('enter name addition [%s_...]: ', ...
    ar.config.optimizers{ar.config.optimizer}), 's');
if(~isempty(name))
    name = sprintf('%s_%s', ar.config.optimizers{ar.config.optimizer}, name);
else
    name = ar.config.optimizers{ar.config.optimizer};
end

ar.fit_hist(j).hist = ar.fit;
ar.fit_hist(j).optimizer = ar.config.optimizer;
ar.fit_hist(j).config = ar.config.optim;
ar.fit_hist(j).name = name;