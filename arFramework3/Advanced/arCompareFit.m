% arCompareFit(indexes, rightalign)
% 
% If two fits are performed, convergence is compared.
% 
%   indexes         
%   rightalign      
function arCompareFit(indexes, rightalign)

global ar

if(~exist('indexes','var') || isempty(indexes))
    indexes = 1:length(ar.fit_hist);
end
if(~exist('rightalign','var'))
    rightalign = false;
end

minchi2 = Inf;
for j=1:length(indexes)
    minchi2 = min([minchi2 ar.fit_hist(indexes(j)).hist.chi2_hist]);
end

minconstr = Inf;
for j=1:length(indexes)
    minconstr = min([minconstr ar.fit_hist(indexes(j)).hist.constr_hist]);
end

minchi2constr = Inf;
for j=1:length(indexes)
    chi2constr = ar.fit_hist(indexes(j)).hist.chi2_hist + ar.fit_hist(indexes(j)).hist.constr_hist;
    minchi2constr = min([minchi2constr chi2constr]);
end

figure(1); clf;
nsub = sum([ar.ndata>0 ar.nconstr>0])+1;
isub = 1;

markerstyle = [];
if(length(indexes)>99)
    markerstyle = 'none';
end

if(ar.ndata>0)
    subplot(1,nsub,isub);
    
    h = [];
    labels = {};
    for j=1:length(indexes)
        C = arLineMarkersAndColors(indexes(j),length(indexes),[],markerstyle,'-');
        qnonnan = ~isnan(ar.fit_hist(indexes(j)).hist.chi2_hist);
        if(~rightalign)
            xs = 1:sum(qnonnan);
        else
            xs = (1:sum(qnonnan)) -1-sum(qnonnan);
        end
        if(sum(qnonnan)>0)
            h(end+1) = semilogy(xs, ar.fit_hist(indexes(j)).hist.chi2_hist(qnonnan) + 1 - minchi2, C{:}); %#ok<AGROW>
            labels{end+1} = ar.fit_hist(indexes(j)).name; %#ok<AGROW>
        end
        hold on
    end
    hold off
    
    if(length(indexes)<21)
        legend(h, strrep(labels, '_', '\_'));
    end
    title('likelihood');
    xlabel('fit iteration');
    grid on
    arSpacedAxisLimits
    
    isub = isub + 1;
end

if(ar.nconstr>0)
    subplot(1,nsub,isub);
    
    for j=1:length(indexes)
        C = arLineMarkersAndColors(indexes(j),length(indexes),[],markerstyle,'-');
        qnonnan = ~isnan(ar.fit_hist(indexes(j)).hist.constr_hist);
        if(~rightalign)
            xs = 1:sum(qnonnan);
        else
            xs = (1:sum(qnonnan)) -1-sum(qnonnan);
        end
        if(sum(qnonnan)>0)
            semilogy(xs, ar.fit_hist(indexes(j)).hist.constr_hist(qnonnan) + 1 - minconstr, C{:});
        end
        hold on
    end
    hold off
    
    title('constraints');
    xlabel('fit iteration');
    grid on
    arSpacedAxisLimits
    
    isub = isub + 1;
end

if(ar.ndata>0 && ar.nconstr>0)  
    subplot(1,nsub,isub);
    
    for j=1:length(indexes)
        chi2constr = ar.fit_hist(indexes(j)).hist.chi2_hist + ar.fit_hist(indexes(j)).hist.constr_hist;
        C = arLineMarkersAndColors(indexes(j),length(indexes),[],markerstyle,'-');
        qnonnan = ~isnan(chi2constr);
        if(~rightalign)
            xs = 1:sum(qnonnan);
        else
            xs = (1:sum(qnonnan)) -1-sum(qnonnan);
        end
        if(sum(qnonnan)>0)
            semilogy(xs, chi2constr(qnonnan) + 1 - minchi2constr, C{:});
        end
        hold on
    end
    hold off
    
    title('likelihood + constraints');
    xlabel('fit iteration');
    grid on
    arSpacedAxisLimits
end

