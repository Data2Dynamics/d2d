%% source DDmore (http://repository.ddmore.foundation/model/DDMODEL00000127#Overview)
% initialize model
arInit
arLoadModel('mod127');
arLoadData('data127');
arCompileAll;

% change prior to regard Kgrw1 as individual
ar.p(7) = 0;
ar.lb(7) = -1;
ar.ub(7) = 2;
ar.ub(2) = 5;

% initial conditions for dep
listAMT = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 5.27, 5.07, 5.27, 5.27, 4.46, 5.48, 5.28, 5.88];
ind = find(~cellfun(@isempty,regexp(ar.pLabel,'p_ID_AMT')));
ar.p(ind) = log10(listAMT(1:length(ind)));
ar.qFit(ind) = 0;

% start values for the optimizer from ddmore
ParamInitPop = [log10(0.0322), log10(111), log10(0.137), log10(6.17), log10(39)];
ParamInitIDsd = [0, 0];
ar.p(1:5) = ParamInitPop;
ar.p(7:8) = ParamInitIDsd;


