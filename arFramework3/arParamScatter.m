function arParamScatter(ncols, nrows)
    global ar
    
    if(isempty(ar.chi2s_sorted))
        arPlotChi2s
    end
    
    sumples = size(ar.pLabel);
    sumples = sumples(2);
    
    if(~exist('ncols', 'var') || isempty(ncols))
        ncols = ceil(sumples^(0.4))+1;
        nrows = ceil(sumples/ncols);
    end
    if(~exist('nrows', 'var') || isempty(nrows))
        nrows = ceil(sumples/ncols);
    end
    somebound = 1;
    min_chi = ar.chi2s_sorted(1);
    chi_step = find(ar.chi2s_sorted < min_chi + somebound);
    ar.chi2s_sorted(chi_step);
    if size(ar.chi2s_sorted(chi_step)) < 5
        warning('there are less than 5 points on the lowest step')
    end
    i = 1;
    j = 1;
    for param = ar.pLabel
        subplot(ncols, nrows, i)
        plot(ar.ps_sorted(1:chi_step(end), i), ar.chi2s_sorted(chi_step), '*')
        title(param)
        xlim([min(ar.ps_sorted(1:chi_step(end), i)) - 1,  max(ar.ps_sorted(1:chi_step(end), i)) + 1])
        ylim([min(ar.chi2s_sorted(chi_step)) - 3,  max(ar.chi2s_sorted(chi_step)) + 3])
        i = i + 1;
    end
end
