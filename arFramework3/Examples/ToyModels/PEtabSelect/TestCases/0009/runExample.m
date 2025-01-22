%% load PEtab Select problem and perform model selection
arPetsSelectModel('petab_select_problem.yaml', d2dFitFunction=@fitFunction);


% -------------------------------------------------------------------------
%% Define custom calibration function
function fitFunction()

global ar

%% Modify d2d configs:
% deactivate Bessel correction
ar.config.useFitErrorCorrection = 0;

% change initial time of model equilibration (0 instead of -1e7)
arClearEvents();
arSteadyState(1, 1, 1, {}, 0);  % fourth argument is initial time

%% perform multi-start fits
arFitLHS(25);

end