%initialize framework, load model and data 
Setup_SIR_EBS

%%
% sinlge fit
arFit;

%% The following code indicates Multi-Start fitting, Profile likelihood calculation and Uncertainty analysis
% % multistart fit
% arFitLHS(100,123);
% % or load best fit parameters
arLoadPars('BestFit_EBS');
% 
% % show parameter values
% arPrint;
% 
% % plot results of multistart fit
% arPlotFits
% 
% % plot observables, internal sates and fluxes
% arQplot('xyv',1);
% arPlot;
% 
% % save results
% arSave('SIR_EBS_Multistart100')
% 
% % likelihood profiles
% arPLEInit
% 
% % Set tolerances
% ar.ple.relchi2stepincrease(5) = 0.01;
% ar.ple.minstepsize(:) = 1e-4;
% 
% % calculate profiles
% ple(1:5,200)
% 
% % plot profiles
% plePlotMulti;
% 
% % plot trajectories along profiles
% arPLETrajectories;
