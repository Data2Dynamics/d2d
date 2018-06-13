%%
close all
Setup


%% kein L1
arLoad(ws_compiled)
global ar pleGlobals

arFit
arFitTillConv

basepath = arSave('NoL1_fitted_log_1e-10');
[pfad,ws_noL1]=fileparts(basepath);
arPlot(true)
arPrint
close all

arLoad(ws_noL1)
pleGlobals.allowbetteroptimum = 1;
arPLEInit

ple([strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)]);
arSave(ws_noL1)

plePlotMulti([],true)
ar.model.qPlotXs(:)=1;
arPlot(true)
ar.model.qPlotXs(:)=0;
arReport

indP = [strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)];
for i=1:length(indP)
    pleExtend(indP(i),50);
    pleExtend(indP(i),-50);
end

%% L1
arLoad(ws_compiled)
global ar pleGlobals

ar = L1penalty(ar,2,1,-10); % log:
arPrint

arFit
arFitTillConv

basepath = arSave('L1_fitted_log_1e-10');
[pfad,ws_fitted]=fileparts(basepath);
arPlot(true)
arPrint
close all

% PLE log
% arLoad(ws_fitted)
pleGlobals.allowbetteroptimum = 1;
arPLEInit

ple([strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)]);
% arSave(ws_fitted)
plePlotMulti([],true)
ar.model.qPlotXs(:)=1;
arPlot(true)
ar.model.qPlotXs(:)=0;
arReport

%% Complexe mit Param am Rand raus:
raus = [10,14,15,16,18,19];
ar.mean(raus) = log10(eps);
ar.p(raus) = log10(eps);
ar.lb(raus) = log10(eps);
ar.qFit(raus) = 0;
arFit
arFitTillConv
arPrint
arSave('NewBounds4_6raus')
arPLEInit
ple(setdiff([strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)],raus));
% print -dpdf PLE_4Nov14_raus
% saveas(gcf,'PLE_4Nov14_raus')

%% Profiles without prior:
plePlotWithoutPrior

%% Verschiedene SDs
% sds = [.2,.5,1,2,5,10];
[ms,sds] = meshgrid(linspace(-10,-4,6),logspace(-2,7,10));
ars = cell(size(sds));

arLoad(ws_fitted);
arIn = ar;

for m=1:size(sds,1)
    for s=1:size(sds,2)
%     arLoad('20140813T115116_L1_PLE_log_1e-7')
        ar = arIn;

        ar = L1penalty(ar,sds(m,s),1,ms(m,s)); % log:

        arFitTillConv
        ars{m,s} = ar;
        
        [ps,Ps] = PlotL1(indP,ms,sds,ars);
    end
end

%%

[ps,Ps] = PlotL1(indP,ms,sds,ars)
saveas(gcf,[ars{1}.config.savepath,'Prior_mean_and_sd2']);


%% r?ckw?rts:
ars2 = ars;
% ars2 = cell(size(ars));
for m=1:size(ars,1)
    for s=(size(ars,2)-1):-1:1
        ar = arIn;
        
        ar = L1penalty(ar,sds(m,s),1,ms(m,s)); % log:
        ar.p = ars{s+1}.p;
        
        arFitTillConv
        ars2{m,s} = ar;
        [ps,Ps] = PlotL1(indP,ms,sds,ars2);
  end
  ars2{m,size(ars,2)} = ars{m,size(ars,2)};
end

%%
[ps,Ps] = PlotL1(indP,ms,sds,ars2)
saveas(gcf,[ars{1}.config.savepath,'/Prior2_mean_and_sd']);


%%

for ip=1:length(indP)
    subplot(3,4,ip)
    surf(log10(sds),ms,Ps{ip});
    axis tight
    title(strrep(ar.pLabel{indP(ip)},'_','\_'))
    set(gcf,'Renderer','zbuffer')
end
%%
    
im = 1;

indP = [strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)];
ps = NaN(size(ars,1),length(indP));
ps2 = NaN(size(ars,1),length(indP));
for i=1:size(ars,1)
    ps(i,:) = ars{i,im}.p(indP);
    ps2(i,:) = ars2{i,im}.p(indP);
end
semilogx(sds(:,im),ps,'.-')
hold on
semilogx(sds(:,im),ps2,'--')
legend(strrep(ar.pLabel(indP),'_','\_'));
ylabel('\theta_{est}')
xlabel('Penalty SD');
title(['Prior mean = ',num2str(ms(im))])

for i=1:size(ps2,2)
    text(sds(end),ps2(end,i),strrep(ar.pLabel{indP(i)},'_','\_'));
end

saveas(gcf,[ar.config.savepath,'/','PenaltyDependency.fig']);
save Setup_ModelSelectionProfile2_L1_1eM10

%% LHS
% 
arFitLHS(100);
arSave('L1_LHS_log_1e-10');
arPlot(true)
arPrint
close all
% 

%%
arLoad('20141104T114117_NewTrdog');
global ar

ar = L1penalty(ar,2,1,-14); % log:

ar.p(16)     = 0; 
ar.qLog10(16)= 0;
ar.qFit(16)  = 0;

ar.p(9) = ar.lb(9);

arFit
arFitTillConv
arSave('without_k_on_u')

arPLEInit
ple([strmatch('k_',ar.pLabel);strmatch('khomo',ar.pLabel)]);
pleExtend
