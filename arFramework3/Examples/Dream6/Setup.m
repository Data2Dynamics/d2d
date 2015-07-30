clear all
%% loading model and a wild-type data-set where everything is observed
disp('arInit')
arInit;

%% Choose on of the three models:

% arLoadModel('model1_dream6');
% arLoadData('all_model1', 1, 'xls',false);      % everything observed

% arLoadModel('model2_dream6');
% arLoadData('all_model2', 1, 'xls',false);        % everything observed

arLoadModel('model3_dream6');
arLoadData('all_model3', 1, 'xls',false);      % everything observed

arCompileAll

%% Loading of the original parameters and simulation of data:
loadGoldStandardParameter(ar.model(1).name) % no priors
arSimuData
arPlot

%% Simulating all single pertubations:
indper = [find(~cellfun(@isempty,regexp(ar.pLabel,'_ko$'))),find(~cellfun(@isempty,regexp(ar.pLabel,'_kd$'))),find(~cellfun(@isempty,regexp(ar.pLabel,'_ic$')))];
ps = ones(length(indper),1)*ar.p;
for i=1:length(indper)
    ps(i,indper(i)) = 1;
end
arPlotMulti(ps)



