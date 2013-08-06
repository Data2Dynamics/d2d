% Initialize and clear workspace of framework
%
% Copyright D2D Development Team
% Contact: Andreas Raue (andreas.raue@fdm.uni-freiburg.de)

arCheck;

global ar

ar = struct([]);
ar(1).stop = 0;
ar.fevals = 0; 

arInitUser;

ar.info.initTime = now;
[ar.info.def_version_code ar.info.c_version_code] = arGetVersion;
fprintf('arFramework3 (.def version %i, c-code version %s)\n', ...
    ar.info.def_version_code, ar.info.c_version_code);
fprintf('Copyright 2013 D2D Development Team. All rights reserved.\n');
fprintf('Contact: Andreas Raue - andreas.raue@fdm.uni-freiburg.de\n\n');

ar.checksum = [];

ar.config.useParallel = true;
ar.config.nParallel = feature('numCores');

% plotting options
ar.config.nFinePoints = 300;
ar.config.savepath = [];
ar.config.par_close_to_bound = 0.1;
ar.config.ssa_min_tau = 1e-3;
ar.config.ssa_runs = 1;

% ar.config.fiterrors:
%   0 = error model (fixed)
%   1 = error model (estimated)
%  -1 = given SD (fixed)
ar.config.fiterrors = 1;
ar.config.fiterrors_correction = 1;
ar.config.fiterrors_correction_warning = false;
ar.config.useFitErrorCorrection = true;

% ar.config.ploterrors:
%   0 = error bands
%   1 = error bars
%  -1 = confidence bands 
ar.config.ploterrors = 0;

ar.config.plotter_sorted = false;

% optimization options
ar.config.useSensis = true; 
ar.config.useJacobian = true;

ar.config.optimizer = 1;
ar.config.optimizers = {'lsqnonlin', 'fmincon', 'empty', 'STRSCNE', 'arNLS'};
ar.config.optim = optimset('lsqnonlin');
ar.config.optim.Display = 'off';
ar.config.optim.TolFun = 0;
ar.config.optim.TolX = 1e-6;
ar.config.optim.MaxIter = 1000;

ar.config.showFitting = 0;

% PPL option structure
ar.ppl.alpha_level = 0.05;
ar.ppl.ndof = 1;
ar.ppl.qLog10 = 1;
ar.ppl.rel_increase = 0.2;
ar.ppl.n = 300;

% CVODES settings
ar.config.atol = 1e-6;
ar.config.rtol = 1e-6;
ar.config.maxsteps = 1000;

ar.config.steady_state_constraint = 1;

% CVODES flags
ar.info.arsimucalc_flags = cell(1,30);
for j = 1:30
    ar.info.arsimucalc_flags{j} = 'not defined';
end
ar.info.arsimucalc_flags{1} = 'malloc UserData';
ar.info.arsimucalc_flags{2} = 'N_VNew_Serial()';
ar.info.arsimucalc_flags{3} = 'CVodeCreate()';
ar.info.arsimucalc_flags{4} = 'CVodeInit()';
ar.info.arsimucalc_flags{5} = 'CVodeSStolerances()';
ar.info.arsimucalc_flags{6} = 'CVodeSetUserData()';
ar.info.arsimucalc_flags{7} = 'CVDense()';
ar.info.arsimucalc_flags{8} = 'CVDlsSetDenseJacFn()';
ar.info.arsimucalc_flags{9} = 'N_VCloneVectorArray_Serial()';
ar.info.arsimucalc_flags{10} = 'CVodeSensInit1()';
ar.info.arsimucalc_flags{11} = 'CVodeSensEEtolerances()';
ar.info.arsimucalc_flags{12} = 'CVodeSetSensErrCon()';
ar.info.arsimucalc_flags{13} = 'CVodeSetSensParams()';
ar.info.arsimucalc_flags{14} = 'CVodeGetSens()';
ar.info.arsimucalc_flags{15} = 'CVodeSetMaxNumSteps';

ar.info.cvodes_flags = cell(1,30);

ar.info.cvodes_flags{1} = 'CV_TOO_MUCH_WORK';
ar.info.cvodes_flags{2} = 'CV_TOO_MUCH_ACC';
ar.info.cvodes_flags{3} = 'CV_ERR_FAILURE';
ar.info.cvodes_flags{4} = 'CV_CONV_FAILURE';

ar.info.cvodes_flags{5} = 'CV_LINIT_FAIL';
ar.info.cvodes_flags{6} = 'CV_LSETUP_FAIL';
ar.info.cvodes_flags{7} = 'CV_LSOLVE_FAIL';
ar.info.cvodes_flags{8} = 'CV_RHSFUNC_FAIL';
ar.info.cvodes_flags{9} = 'CV_FIRST_RHSFUNC_ERR';
ar.info.cvodes_flags{10} = 'CV_REPTD_RHSFUNC_ERR';
ar.info.cvodes_flags{11} = 'CV_UNREC_RHSFUNC_ERR';
ar.info.cvodes_flags{12} = 'CV_RTFUNC_FAIL';

ar.info.cvodes_flags{20} = 'CV_MEM_FAIL';
ar.info.cvodes_flags{21} = 'CV_MEM_NULL';
ar.info.cvodes_flags{22} = 'CV_ILL_INPUT';
ar.info.cvodes_flags{23} = 'CV_NO_MALLOC';
ar.info.cvodes_flags{24} = 'CV_BAD_K';
ar.info.cvodes_flags{25} = 'CV_BAD_T';
ar.info.cvodes_flags{26} = 'CV_BAD_DKY';
ar.info.cvodes_flags{27} = 'CV_TOO_CLOSE';



