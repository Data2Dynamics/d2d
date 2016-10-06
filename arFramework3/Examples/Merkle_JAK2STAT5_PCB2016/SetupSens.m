%% Load models & data

arInit;

arLoadModel('jak2_stat5_cfue_sens');
arLoadModel('jak2_stat5_h838_l1_final_sens');

arCompileAll;

arLoadPars('sensitivityanalysis')
ar.lb = 10.^ar.lb;
ar.ub = 10.^ar.ub;

% arSetPars(pLabel, p, qFit, qLog10, lb, ub, type, meanp, stdp)
arSetPars('ActD',0,0,0,0,0)
arSetPars('SOCS3oe',0,0,0,0,0)
arSetPars('overexp',0,0,0,0,0)

arSetPars('epo_level',-6,0,1,-8,-4)

arFindInputs

%% Sensitivity analysis
ar.p(ar.qLog10==1) = 10.^ar.p(ar.qLog10==1);
ar.qLog10(ar.qLog10==1) = 0;

ip1 = [3:19 21 23:28 33];
ip2 = [3:14 17:21 23 26:31 37];

pNames1 = ar.model(1).condition(1).p(ip1);
pNames2 = ar.model(2).condition(1).p(ip2);

for i = 1:length(pNames1)
    if ~strcmp(pNames1{i},strrep(pNames2{i},'_H838',''))
        error('Parameter do not match!!')
    end
end

npSTAT5 = strcmp(ar.model(1).xNames,'npSTAT5_auc');
[a,b] = min(abs(ar.model(1).condition(1).tFine-60));
[c,d] = min(abs(ar.model(2).condition(1).tFine-60));

arSimu(1,1) % Sens + Fine

sens_cfue = ar.model(1).condition(1).pNum(ip1)' / ar.model(1).condition(1).xFineSimu(b,npSTAT5) .* squeeze(ar.model(1).condition(1).sxFineSimu(b,npSTAT5,ip1));
sens_h838 = ar.model(2).condition(1).pNum(ip2)' / ar.model(2).condition(1).xFineSimu(d,npSTAT5) .* squeeze(ar.model(2).condition(1).sxFineSimu(d,npSTAT5,ip2));

figure
[sens_sort,c] = sort(sens_h838);
barh(sens_sort)
set(gca,'YTick',1:length(ip1),'YTickLabel',pNames1(c),'YLim',[0 length(ip1)+1],'XLim',[-1 1])

figure
barh([sens_h838(c) sens_cfue(c)])
set(gca,'YTick',1:length(ip1),'YTickLabel',pNames1(c),'YLim',[0 length(ip1)+1],'XLim',[-1 1])

y = linspace(-1,1,201);
colmap = nan(length(y),3);

for iy = 1:length(y)
    if iy < length(y)/2
        sat = 1-iy/(length(y)/2);
        hue = 2/3;
        val = 1;
    elseif iy == length(y)/2-.5
        sat = 0;
        hue = 0;
        val = 1;
    else
        sat = (iy-length(y)/2)/length(y);
        hue = 0;
        val = 1;
    end
    colmap(iy,:) = hsv2rgb([hue sat val]);
end

figure
for ix = 1:length(sens_h838)
    patch([ix ix+1 ix+1 ix],[0 0 1 1],ones(1,4),'FaceColor',colmap(max([1 round((sens_h838(ix)+1)/2*201)]),:),'EdgeColor','none')
end
set(gca,'XLim',[1 length(sens_h838)+1])
set(gca,'XTick',(1:length(ip1))+.5,'XTickLabel',pNames1)