%% Load results
% Skip the 5 already existing workspaces. Only workspaces with filenames
% containing 'result' should be loaded. Skip failed L1 scans.
existing = 5;
failed = [71 97 147 153 162 186 201 299 328 337 457 476];
oks = nSimu + existing + (setdiff(1:nSimu,failed));

ws = {};
for i = 1:length(oks)
    arLoad(oks(i))
    ws{end+1} = ar;
end

%% Plot ROC curve
fpr = nan(length(ar.L1linv)+1,length(oks));
tpr = nan(length(ar.L1linv)+1,length(oks));
fprfinal = nan(1,length(oks));
tprfinal = nan(1,length(oks));
for i = 1:length(oks)
    ar = ws{i};
    ar.L1ps = [ar.L1ps; zeros(1,size(ar.L1ps,2))];
    
    psTrue = repmat(ar.pTrue(relto),length(ar.L1linv)+1,1);
    tn = ar.L1ps(:,relto) == 0 & psTrue == 0;
    tp = ar.L1ps(:,relto) ~= 0 & psTrue ~= 0;
    fn = ar.L1ps(:,relto) == 0 & psTrue ~= 0;
    fp = ar.L1ps(:,relto) ~= 0 & psTrue == 0;
    
    fpr(:,i) = sum(fp,2)./(sum(fp,2) + sum(tn,2));
    tpr(:,i) = sum(tp,2)./(sum(tp,2) + sum(fn,2));
    
    fprfinal(i) = fpr(ar.L1final_ind,i);
    tprfinal(i) = tpr(ar.L1final_ind,i);
end

%%
fprfine = 0:0.01:1;
tprfine = nan(length(fprfine),size(fpr,2));

for i = 1:length(oks)
    [fpr(:,i),b] = sort(fpr(:,i));
    tpr(:,i) = tpr(b,i);
    [unif,unit] = unique(fpr(:,i));
    for j = 1:length(unif)
        unit(j) = mean(tpr(fpr(:,i) == unif(j),i));
    end
    
    tprfine(:,i) = interp1(unif,unit,fprfine,'linear','extrap');
end

fpruni = fprfine';
tpruni = mean(tprfine,2);
tprunistd = std(tprfine,[],2);

%%
chi2 = nan(size(oks));
for i = 1:length(oks)
    chi2(i) = sum((tprfine(:,i)-tpruni).^2);
end
[chi2_sorted,sorti] = sort(chi2);
% The ROC curves closest to the mean:
% --> [287   344   347     8   278   337   391    12   263   156 ...]

%%

figure
plot(fpruni,tpruni,'k')
hold on
patch([fpruni; fpruni(end:-1:1)],[tpruni-tprunistd; tpruni(end:-1:1)+tprunistd(end:-1:1)],ones(2*length(fpruni),1),'FaceColor','k','FaceAlpha',.2)
plot(fprfinal(sorti(1)),tprfinal(sorti(1)),'ro')
plot(fprfine,tprfine(:,sorti(1)),'r')
xlabel('FPR')
ylabel('TPR')
title('ROC')
set(gca,'XLim',[0 1],'YLim',[0 1])
plot([0 1],[0 1],'k--')

%%

figure
plot(fpruni,tpruni,'k')
hold on
patch([fpruni; fpruni(end:-1:1)],[tpruni-tprunistd; tpruni(end:-1:1)+tprunistd(end:-1:1)],ones(2*length(fpruni),1),'FaceColor','k','FaceAlpha',.2)
plot(fprfinal,tprfinal,'o')
xlabel('FPR')
ylabel('TPR')
title('ROC')
set(gca,'XLim',[0 1],'YLim',[0 1])
plot([0 1],[0 1],'k--')

%%
pTrues = nan(length(ws),length(ar.p));
pSels = pTrues;
for i = 1:length(ws)
    pTrues(i,:) = ws{i}.pTrue;
    pSels(i,:) = ws{i}.p;
end

tnall = pSels == 0 & pTrues == 0;
tpall = pSels ~= 0 & pTrues ~= 0;
fnall = pSels == 0 & pTrues ~= 0;
fpall = pSels ~= 0 & pTrues == 0;

strs = {'p_degradation_rate','pro','rbs','Kd','_h'};
fprCat = nan(1,length(strs));
tprCat = fprCat;
nCat = fprCat;
accCat = fprCat;
nTest = [500 3000 3000 4000 4000];

fcsFine = [2 4 5 10];
fcsFine = [1./fcsFine fcsFine];
fcsFine = log10(fcsFine);

fprCatFine = nan(length(strs),length(fcsFine));
tprCatFine = fprCatFine;
nCatFine = fprCatFine;
for i = 1:length(strs)
    parCat = ~cellfun(@isempty,strfind(ar.pLabel,strs{i}));
    fprCat(i) = sum(sum(fpall(:,relto(parCat(relto)))))/(sum(sum(fpall(:,relto(parCat(relto))))) + sum(sum(tnall(:,relto(parCat(relto))))));
    tprCat(i) = sum(sum(tpall(:,relto(parCat(relto)))))/(sum(sum(tpall(:,relto(parCat(relto))))) + sum(sum(fnall(:,relto(parCat(relto))))));
    nCat(i) = sum(sum(pTrues(:,relto(parCat(relto)))~=0));
    accCat(i) = (sum(sum(tpall(:,relto(parCat(relto))))) + sum(sum(tnall(:,relto(parCat(relto)))))) / nTest(i);
    
    parCat = repmat(parCat,size(pTrues,1),1);
    reltos = parCat*0;
    reltos(:,relto) = 1;
    for j = 1:length(fcsFine)
        myp = (pTrues == fcsFine(j) | pTrues == 0) & reltos & parCat;
        fprCatFine(i,j) = sum(sum(fpall(myp)))/(sum(sum(fpall(myp))) + sum(sum(tnall(myp))));
        tprCatFine(i,j) = sum(sum(tpall(myp)))/(sum(sum(tpall(myp))) + sum(sum(fnall(myp))));
        nCatFine(i,j) = sum(sum(pTrues(myp)~=0));
    end
end

fprAll = sum(sum(fpall(:,relto)))/(sum(sum(fpall(:,relto))) + sum(sum(tnall(:,relto))));
tprAll = sum(sum(tpall(:,relto)))/(sum(sum(tpall(:,relto))) + sum(sum(fnall(:,relto))));
nAll = sum(sum(pTrues(:,relto)~=0));
accAll = (sum(sum(tpall(:,relto))) + sum(sum(tnall(:,relto))))/(length(ws)*length(relto));


%%
relto = arPrint('relto');
parCat = ~cellfun(@isempty,strfind(ar.pLabel,strs{5}));
relto = relto(parCat(relto));

fcs = unique(pTrues(:,relto));

pest = pSels(:,relto);
ptru = pTrues(:,relto);
ptru = ptru(:);
pest = pest(:);
neqzero = pest~=0;
ptru = ptru(neqzero);
pest = pest(neqzero);

closest = nan(size(pest));
signs = nan(size(pest));
for i = 1:length(pest)
    [a,b] = min(abs(pest(i)-fcs));
    closest(i) = ptru(i) == fcs(b);
    signs(i) = abs(sign(pest(i))-sign(ptru(i))) < 2;
end

accCloseH = sum(closest)/length(closest);
accSignH = sum(signs)/length(signs);

relto2 = arPrint('relto');
relto = setdiff(relto2,relto);

fcs = unique(pTrues(:,relto));

pest = pSels(:,relto);
ptru = pTrues(:,relto);
ptru = ptru(:);
pest = pest(:);
neqzero = pest~=0;
ptru = ptru(neqzero);
pest = pest(neqzero);

sumH = length(closest);
closest = [closest; nan(size(pest))];
signs = [signs; nan(size(pest))];
for i = 1:length(pest)
    [a,b] = min(abs(pest(i)-fcs));
    closest(i+sumH) = ptru(i) == fcs(b);
    signs(i+sumH) = abs(sign(pest(i))-sign(ptru(i))) < 2;
end

accClose = sum(closest)/length(closest);
accSign = sum(signs)/length(signs);

%%
arLoad('1_result_294')
relto = arPrint('relto');

unsup = sum(ar.p(relto)'==0 & ar.pTrue(relto)'==0 | ar.p(relto)'~=0 & ar.pTrue(relto)'~=0)/length(relto);
sup = (sum(ar.p(relto)'==0 & ar.pTrue(relto)'==0 | ar.p(relto)'~=0 & ar.pTrue(relto)'~=0)+2)/length(relto);
hypoth = (sum(ar.p(relto)'==0 & ar.pTrue(relto)'==0 | ar.p(relto)'~=0 & ar.pTrue(relto)'~=0)+4)/length(relto);