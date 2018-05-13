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

clear all;
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
arSimu(true,true,true);

%% Parameter settings
% Note: The parameter values and bounds are set to reproduce the
% simulations in the original publication.

% Parameters with linear influence
arSetPars('p11',-0.18931146441630517 ,[],[], -5, 5)  % PKD phosphorylation influenced by Caramide inflow
arSetPars('p12',-1.2951693058037776  ,[],[], -5, 5)  % PKD phosphorylation
arSetPars('p13',-0.68528208807529756 ,[],[], -5, 5)  % PKD dephosphorylation
arSetPars('p21',0.91322139750063747  ,[],[], -5, 5)  % PI4K3B dephosphorylation
arSetPars('p22',0.26999507129007366  ,[],[], -5, 5)  % PI4K3Ba phosphorylation by PKD
arSetPars('p31',2.8038329601560896   ,[],[], -5, 5)  % CERT recrution to the TGN membrane by aPI4K3b
arSetPars('p32',0.67759061446638846  ,[],[], -5, 5)  % CERT dephosphorylation (activation) at ER
arSetPars('p33',4.1092688852009367   ,[],[], -5, 5)  % CERT phosphorylation (deactivation) by PKDDAGa and detachment from TGN

% Degradation rates
arSetPars('a11',-0.002926956967776128,[],[], -4, 2)  % PKD basal degradation
arSetPars('a12',0.96790906144441557  ,[],[], -4, 2)  % PKDa basal degradation
arSetPars('a21',0.11866927487115031  ,[],[], -4, 2)  % PI4K3B degradation
arSetPars('a22',-2.5206963464088994  ,[],[], -4, 2)  % PI4K3Ba degradation
arSetPars('a31',-3.6407643715083622  ,[],[], -4, 2)  % CERT degradation rates
arSetPars('a32',-3.1066933752303902  ,[],[], -4, 2)  % CERT degradation rates
arSetPars('a33',-3.4914558663191695  ,[],[], -4, 2)  % CERT degradation rates

% Synthesis rates
arSetPars('s12',5.6858407186167481   ,[],[],  2, 7)  % PKD basal production rate
arSetPars('s21',6.3200257512233726   ,[],[],  2, 7)  % PI4K3B basal production rate
arSetPars('s31',2.1432339301626024   ,[],[],  2, 7)  % CERT basal production rate

% Input scalings
arSetPars('pu2',0,0)                                   % ectopical expression of PKD strength
arSetPars('pu3',7.9964407645546816   ,[],[], -1, 8)    % ectopical expression of PI4K3B strength
arSetPars('pu4',7.1804513037054267   ,[],[], -1, 8)    % ectopical expression of CERT strength
arSetPars('pu5',1.1570329047724215   ,[],[], -1, 8)    % PdBu aktivation strength
arSetPars('pu6',1.3021265927148578   ,[],[], -1, 8)    % kbNB142-70 inhibition strength 

% Threshold parameters
arSetPars('m11',7.3694401379449079   ,[],[],  2,10)  % PKD activation via Caramide inflow saturation threshold
arSetPars('m22',3.8568422695340612   ,[],[],  2,10)  % PI4K3Ba phosphorylation saturation constant
arSetPars('m31',9.6534422825120316   ,[],[],  2,10)  % CERT phosphorylation saturation
arSetPars('m33',6.8860692749081878   ,[],[],  2,10)  % CERT dephosphorylation saturation

% Scaling factors and standard deviations which have been introduced to 
% avoid division of state variables using by Weber et al.
% (Note: The values for the scaling are set to test the implementation.)
arSetPars('scale_yPKDpN0'     ,log10(1/2.4931e+03))
arSetPars('scale_yPKDpN24'    ,log10(1/2.8950e+03))
arSetPars('scale_yPKDpN25'    ,log10(1/2.6302e+04))
arSetPars('scale_yPI4K3BpRN24',log10(1/0.0612)    )
arSetPars('scale_yCERTpRN24'  ,log10(1/0.0062)    )

% (Note: The values of the standard deviation of the absolute measurements
% are set as they are condensed in a single replicate.)
arSetPars('std_yPKDt'       ,log10(  1910527.835),0,[],-2,10)
arSetPars('std_yPI4K3Bt'    ,log10(   887676.131),0,[],-2,10)
arSetPars('std_yCERTt'      ,log10(244727220.419),0,[],-2,10)
arSetPars('std_yPKDpN0'     ,[]   ,[],[],-2,2)
arSetPars('std_yPKDpN24'    ,[]   ,[],[],-2,2)
arSetPars('std_yPKDpN25'    ,[]   ,[],[],-2,2)
arSetPars('std_yPI4K3BpRN24',[]   ,[],[],-2,2)
arSetPars('std_yCERTpRN24'  ,[]   ,[],[],-2,2)

%% Parameter optimization
arFitLHS(200)
arPlotChi2s

%% Visualization
arSimu(true,true,true);
arPlot;

