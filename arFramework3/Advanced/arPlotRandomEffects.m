% plot random effect distributions

function arPlotRandomEffects

global ar

if isfield(ar,'random')
    
    figure(1);
    
    [nrows, ncols] = arNtoColsAndRows(length(ar.random));
    
    for j = 1:length(ar.random)
        subplot(nrows, ncols, j);
        
        hist(ar.p(ar.random{j}), 50);
        
        label_is = find(ar.random{j});
        title(strrep(ar.pLabel(label_is(1)),'_','\_'));
    end
else
    warning('No random effects defined.')
end

    