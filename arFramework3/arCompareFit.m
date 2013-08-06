function arCompareFit(indexes)

global ar

if(~exist('indexes','var'))
    indexes = 1:length(ar.fit_hist);
end

minchi2 = Inf;
for j=1:length(indexes)
    minchi2 = min([minchi2 ar.fit_hist(j).hist.chi2_hist]);
end

minconstr = Inf;
for j=1:length(indexes)
    minconstr = min([minconstr ar.fit_hist(j).hist.constr_hist]);
end

minchi2constr = Inf;
for j=1:length(indexes)
    chi2constr = ar.fit_hist(j).hist.chi2_hist + ar.fit_hist(j).hist.constr_hist;
    minchi2constr = min([minchi2constr chi2constr]);
end

figure(1); clf;
nsub = sum([ar.ndata>0 ar.nconstr>0])+1;
isub = 1;

if(ar.ndata>0)
    subplot(1,nsub,isub);
    
    h = nan(1,length(indexes));
    labels = {};
    for j=1:length(indexes)
        C = arLineMarkersAndColors(j,length(indexes),[],[],'-');
        qnonnan = ~isnan(ar.fit_hist(j).hist.chi2_hist);
        h(j) = semilogy(ar.fit_hist(j).hist.chi2_hist(qnonnan) + 1 - minchi2, C{:});
        labels{j} = ar.fit_hist(j).name; %#ok<AGROW>
        hold on
    end
    hold off
    
    legend(h, strrep(labels, '_', '\_'));
    title('likelihood');
    xlabel('fit iteration');
    grid on
    
    isub = isub + 1;
end

if(ar.nconstr>0)
    subplot(1,nsub,isub);
    
    for j=1:length(indexes)
        C = arLineMarkersAndColors(j,length(indexes),[],[],'-');
        qnonnan = ~isnan(ar.fit_hist(j).hist.constr_hist);
        semilogy(ar.fit_hist(j).hist.constr_hist(qnonnan) + 1 - minconstr, C{:});
        hold on
    end
    hold off
    
    title('constraints');
    xlabel('fit iteration');
    grid on
    
    isub = isub + 1;
end

if(ar.ndata>0 && ar.nconstr>0)  
    subplot(1,nsub,isub);
    
    for j=1:length(indexes)
        chi2constr = ar.fit_hist(j).hist.chi2_hist + ar.fit_hist(j).hist.constr_hist;
        C = arLineMarkersAndColors(j,length(indexes),[],[],'-');
        qnonnan = ~isnan(chi2constr);
        semilogy(chi2constr(qnonnan) + 1 - minchi2constr, C{:});
        hold on
    end
    hold off
    
    title('likelihood + constraints');
    xlabel('fit iteration');
    grid on
end