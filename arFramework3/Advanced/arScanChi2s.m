function arScanChi2s(jk, N, doYs)

global ar

if(length(jk)>1)
    error('length(jk) != 1');
end
if(~exist('N','var'))
    N = 100;
end

if(~exist('doYs','var'))
    doYs = false;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

pReset = ar.p;
bestp = ar.p(jk);
[bestchi2s, labels] = arPrintChi2(doYs);

arWaitbar(0);

ps = linspace(ar.lb(jk), ar.ub(jk), N);
chi2s = nan(N,length(bestchi2s));    

for j=1:length(ps);
    arWaitbar(j,N);
    ar.p(jk) = ps(j);
        try
            arCalcMerit(false);
            chi2s(j,:) = arPrintChi2(doYs) - bestchi2s;
        catch error_id
            fprintf('%s for %g: %s\n', ar.pLabel{jk}, ps(j), error_id.message);
        end 
end
arWaitbar(-1);

ar.p = pReset;
arCalcMerit(false);

stds = std(chi2s,0,1);
[~, indexes] = sort(stds,2,'descend');
for j=1:length(stds)
    fprintf('%g\t %s\n', stds(indexes(j)), labels{indexes(j)});
end

nmax = min([length(stds) 8]);

figure(1)
plot(ps, chi2s(:,indexes(1:nmax)));
legend(strrep(labels(indexes(1:nmax)),'_','\_'))
hold on
if(size(chi2s,2) > nmax)
    plot(ps, chi2s(:,indexes((nmax+1):end)), '--');
end
plot([bestp bestp], ylim, 'k--');
hold off

if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1) || ...
        (ar.config.useFitErrorMatrix ==1 && sum(sum(ar.config.fiterrors_matrix==1))>0) )
    ylabel('-2*log(L) increase');
else
    ylabel('chi^2 increase');
end
xlabel(strrep(ar.pLabel{jk},'_','\_'))




