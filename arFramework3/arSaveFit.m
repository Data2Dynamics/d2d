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