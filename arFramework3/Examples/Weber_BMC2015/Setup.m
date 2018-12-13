% This routine runs the parameter estimation for the model of PKD and CERT
% interactions at the trans-Golgi network of mammalian cells introduced in
% the publication
%
%   Weber, P., Hornjik. M., Olayioye, M. A., Hausser, A., and Radde, N. E.
%   (2015) A computational model of PKD and CERT interactions at the
%   trans-Golgi network of mammalian cells. BMC Systems Biology, 20159:9.
%   doi: 10.1186/s12918-015-0147-1
%
% The original publication contains observables which are normalized to a
% specific time point, e.g. yPKDpN24(t) = PKDDAGa(t)/PKDDAGa(t=24h), which
% cannot be handled in D2D. To address this problem, we introduced
% introduce scaling constants for these variables and the set the
% measurement at the respective reference point to one. The resulting
% estimation problem is similar to the original one but not identical.

close all;
clc;

%% Compile model
arInit
arLoadModel('model_A');
arLoadData('experiment_0',1); % Condition introduced for steady state equlibration
arLoadData('experiment_1',1);
arLoadData('experiment_2',1);
arCompileAll;

%% Set pre-equilibration
arSteadyState(1,arFindCondition(ar,'experiment_0'),1:length(ar.model.condition));

%% Parameter settings
% Note: The parameter values and bounds are set to reproduce the
% simulations in the original publication.

% Scaling factors and standard deviations which have been introduced to 
% avoid division of state variables using by Weber et al.
% The values of the standard deviation of the absolute measurements
% are set as they are condensed in a single replicate.

arLoadPars('bestFit')
arFindInputs;
%% Visualization
arSimu(true,true,true);
arPlot;

