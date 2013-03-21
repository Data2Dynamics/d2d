% chi^2 test for model to data compliance

function arChi2Test
global ar

ar.pval = 1-chi2cdf(ar.chi2, ar.ndata - sum(ar.qFit==1));
ar.pval_reduced = 1-chi2cdf(ar.chi2, ar.ndata);

if(ar.pval>=ar.config.alpha_level)
    fprintf('pval(%i) = %f, model is compliant with data for %5.2f%% sign. level\n', ...
        ar.ndata - sum(ar.qFit==1), ar.pval, (1-ar.config.alpha_level)*100);
else
    fprintf('pval(%i) = %f, model is NOT compliant with data for %5.2f%% sign. level\n', ...
        ar.ndata - sum(ar.qFit==1), ar.pval, (1-ar.config.alpha_level)*100);
end

if(ar.pval_reduced>=ar.config.alpha_level)
    fprintf('pval(%i) = %f, model is compliant with data for %5.2f%% sign. level\n', ...
        ar.ndata, ar.pval_reduced, (1-ar.config.alpha_level)*100);
else
    fprintf('pval(%i) = %f, model is NOT compliant with data for %5.2f%% sign. level\n', ...
        ar.ndata, ar.pval_reduced, (1-ar.config.alpha_level)*100);
end