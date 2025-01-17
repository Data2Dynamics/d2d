% load PEtab Select problem and perform model selection
arPetsSelectModel('petab_select_problem.yaml', d2dFitFunction=@fitFunction);

function fitFunction()

global ar

% deactivate Bessel correction
ar.config.useFitErrorCorrection = 0;

% % use settings similar to pyPESTO
% ar.config.rtol = 1e-16;
% ar.config.atol = 1e-12;
% ar.config.maxsteps = 1e6;
% 
% % stricter equilibration settings
% ar.config.eq_rtol = 1e-12;
% ar.config.eq_tol = 1e-12;
% ar.config.init_eq_step = 1e3;
% 
% % stricter optimization settings
% ar.config.optim.MaxIter = 1e6;
% ar.config.optim.TolX = 1e-16;

% change steady state calculation start time (0 instead of -1e7)
arClearEvents();
arSteadyState(1, 1, 1, {}, 0);

% perform fits
arFit();
arPetsFitWithSmartInitials();
arFitLHS(10);

end