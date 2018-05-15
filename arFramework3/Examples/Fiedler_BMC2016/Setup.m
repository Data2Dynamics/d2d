clear all;
close all;
clc;

%% Compile model
arInit
arLoadModel('model_RafMekErk');
arLoadData('Ctrl',1);
arLoadData('Sorafenib_5_muM',1);
arLoadData('UO126_30_muM',1);
arCompileAll;

%% Set pre-equilibration
arSteadyState(1,arFindCondition(ar,'Ctrl'),1);

%% Parameter settings
arLoadPars('bestFit')

%% Visualization
arSimu(true,true,true);
arPrint
figure;

maxMEK = max([max(ar.model.data(1).yFineSimu(:,1)/arGetPars('s_pMek_20140430_gel1',0)),...
              max(ar.model.data(2).yFineSimu(:,1)/arGetPars('s_pMek_20140505_gel1',0)),...
              max(ar.model.data(3).yFineSimu(:,1)/arGetPars('s_pMek_20140505_gel1',0))]);
maxERK = max([max(ar.model.data(1).yFineSimu(:,2)/arGetPars('s_pErk_20140430_gel1',0)),...
              max(ar.model.data(2).yFineSimu(:,2)/arGetPars('s_pErk_20140505_gel1',0)),...
              max(ar.model.data(3).yFineSimu(:,2)/arGetPars('s_pErk_20140505_gel1',0))]);

% control - pMEK
subplot(2,3,1)
fill([ar.model.data(1).tFine;ar.model.data(1).tFine(end:-1:1)],...
     [ar.model.data(1).yFineSimu(:       ,1)+ar.model.data(1).ystdFineSimu(:       ,1);...
      ar.model.data(1).yFineSimu(end:-1:1,1)-ar.model.data(1).ystdFineSimu(end:-1:1,1)]/arGetPars('s_pMek_20140430_gel1',0)/maxMEK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
h(1) = plot(ar.model.data(1).tFine,ar.model.data(1).yFineSimu(:,1)/arGetPars('s_pMek_20140430_gel1',0)/maxMEK,'k-');
h(2) = plot(ar.model.data(1).tExp,ar.model.data(1).yExp(:,1)/arGetPars('s_pMek_20140430_gel1',0)/maxMEK,'ko'); hold on; 
h(3) = plot(ar.model.data(1).tExp,ar.model.data(1).yExp(:,3)/arGetPars('s_pMek_20140505_gel1',0)/maxMEK,'k+');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pMek [UI]')

% 5 muM Sorafenib - pMEK
subplot(2,3,2)
fill([ar.model.data(2).tFine;ar.model.data(2).tFine(end:-1:1)],...
     [ar.model.data(2).yFineSimu(:       ,1)+ar.model.data(2).ystdFineSimu(:       ,1);...
      ar.model.data(2).yFineSimu(end:-1:1,1)-ar.model.data(2).ystdFineSimu(end:-1:1,1)]/arGetPars('s_pMek_20140505_gel1',0)/maxMEK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
plot(ar.model.data(2).tFine,ar.model.data(2).yFineSimu(:,1)/arGetPars('s_pMek_20140505_gel1',0)/maxMEK,'k-');
plot(ar.model.data(2).tExp,ar.model.data(2).yExp(:,1)/arGetPars('s_pMek_20140430_gel2',0)/maxMEK,'k>');
plot(ar.model.data(2).tExp,ar.model.data(2).yExp(:,3)/arGetPars('s_pMek_20140505_gel2',0)/maxMEK,'kx');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pMek [UI]')

% 30 muM UO126 - pMEK
subplot(2,3,3)
fill([ar.model.data(3).tFine;ar.model.data(3).tFine(end:-1:1)],...
     [ar.model.data(3).yFineSimu(:       ,1)+ar.model.data(3).ystdFineSimu(:       ,1);...
      ar.model.data(3).yFineSimu(end:-1:1,1)-ar.model.data(3).ystdFineSimu(end:-1:1,1)]/arGetPars('s_pMek_20140505_gel1',0)/maxMEK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
plot(ar.model.data(3).tFine,ar.model.data(3).yFineSimu(:,1)/arGetPars('s_pMek_20140505_gel1',0)/maxMEK,'k-');
plot(ar.model.data(3).tExp,ar.model.data(3).yExp(:,1)/arGetPars('s_pMek_20140430_gel2',0)/maxMEK,'k>');
plot(ar.model.data(3).tExp,ar.model.data(3).yExp(:,3)/arGetPars('s_pMek_20140505_gel2',0)/maxMEK,'kx');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pMek [UI]')

% control - pERK
subplot(2,3,4)
fill([ar.model.data(1).tFine;ar.model.data(1).tFine(end:-1:1)],...
     [ar.model.data(1).yFineSimu(:       ,2)+ar.model.data(1).ystdFineSimu(:       ,2);...
      ar.model.data(1).yFineSimu(end:-1:1,2)-ar.model.data(1).ystdFineSimu(end:-1:1,2)]/arGetPars('s_pErk_20140430_gel1',0)/maxERK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
h(1) = plot(ar.model.data(1).tFine,ar.model.data(1).yFineSimu(:,2)/arGetPars('s_pErk_20140430_gel1',0)/maxERK,'k-');
h(2) = plot(ar.model.data(1).tExp,ar.model.data(1).yExp(:,2)/arGetPars('s_pErk_20140430_gel1',0)/maxERK,'ko'); hold on; 
h(3) = plot(ar.model.data(1).tExp,ar.model.data(1).yExp(:,4)/arGetPars('s_pErk_20140505_gel1',0)/maxERK,'k+');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pErk [UI]')

% 5 muM Sorafenib - pERK
subplot(2,3,5)
fill([ar.model.data(2).tFine;ar.model.data(2).tFine(end:-1:1)],...
     [ar.model.data(2).yFineSimu(:       ,2)+ar.model.data(2).ystdFineSimu(:       ,2);...
      ar.model.data(2).yFineSimu(end:-1:1,2)-ar.model.data(2).ystdFineSimu(end:-1:1,2)]/arGetPars('s_pErk_20140430_gel2',0)/maxERK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
plot(ar.model.data(2).tFine,ar.model.data(2).yFineSimu(:,2)/arGetPars('s_pErk_20140430_gel2',0)/maxERK,'k-');
plot(ar.model.data(2).tExp,ar.model.data(2).yExp(:,2)/arGetPars('s_pErk_20140430_gel2',0)/maxERK,'k>');
plot(ar.model.data(2).tExp,ar.model.data(2).yExp(:,4)/arGetPars('s_pErk_20140505_gel2',0)/maxERK,'kx');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pErk [UI]')

% 30 muM UO126 - pERK
subplot(2,3,6)
fill([ar.model.data(3).tFine;ar.model.data(3).tFine(end:-1:1)],...
     [ar.model.data(3).yFineSimu(:       ,2)+ar.model.data(3).ystdFineSimu(:       ,2);...
      ar.model.data(3).yFineSimu(end:-1:1,2)-ar.model.data(3).ystdFineSimu(end:-1:1,2)]/arGetPars('s_pErk_20140430_gel2',0)/maxERK,...
      'k','facecolor',0.9*[1,1,1],'edgecolor',0.9*[1,1,1]); hold on; 
plot(ar.model.data(3).tFine,ar.model.data(3).yFineSimu(:,2)/arGetPars('s_pErk_20140430_gel2',0)/maxERK,'k-');
plot(ar.model.data(3).tExp,ar.model.data(3).yExp(:,2)/arGetPars('s_pErk_20140430_gel2',0)/maxERK,'k>');
plot(ar.model.data(3).tExp,ar.model.data(3).yExp(:,4)/arGetPars('s_pErk_20140505_gel2',0)/maxERK,'kx');
axis([0,10,-0.25,1.25])
xlabel('time [h]')
ylabel('relative pErk [UI]')
