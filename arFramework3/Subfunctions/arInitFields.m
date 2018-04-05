% This function checks whether all the required fields in the ar structure
% are present. If some are missing, which occurs when the ar struct was 
% saved with a previous version of D2D, they will be set to the default 
% value.

function ar = arInitFields(ar)
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!
    % !!  NOTE: Every time you add or remove a field, increment this value by one.
    % !! 
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arFormatVersion = 4;
    
    % Without arguments, just return the version number
    if ( nargin < 1 )
        ar = arFormatVersion;
        return;
    end

    % Add substructures if they don't exist
    ar              = checkForField(ar, 'config');
    ar              = checkForField(ar, 'info');
    ar              = checkForField(ar, 'ppl');
    ar              = checkForField(ar, 'ple');
    
    % Config options
    defaults = { ...
        {'fastEquilibration',           false}, ...                     % Faster equilibration (BETA). Set to 1 before compiling if you want to enable rootfinding from within C.
        {'turboSplines',                false}, ...                     % Faster splines (BETA).
        {'turboSSSensi',                false}, ...                     % Faster equilibration (BETA). Toggle with arFastSensis. DO NOT TOGGLE BY HAND.
        {'sensitivitySubset',           0}, ...                         % Only compute subset of sensitivities when certain qFit's are 0 (BETA)
        {'lightSave',                   false}, ...                     % When calling arSave, only save parameter sets by default (useful for big models)
        {'saveMexFile',                true}, ...                      % When calling arSave, the arSimuCalc mex file is copied to the savefolder
        {'skipSim',                     false}, ...                     % Disable simulation (used for fitting steady state models)
        {'barhack',                     false}, ...                     % Display data with only a single time point as bar
        {'networkgraph',                false}, ...                     % Enable network graph plotting features (disabled by default since they slow compilation)
        {'debugExp',                    false}, ...                     % Show crosses with simulation values used for residual calculation in plots
        {'useCache',                    0}, ...                         % Use caching system (0 = no, 1 = strict (also check tExp, tFine), 2 = sloppy)
        {'checkForNegFluxes',           true}, ...
        {'useParallel',                 true}, ...                      % Parallelization
        {'nCore',                       feature('numCores')}, ...       %   number of available cores
        {'nParallel',                   2*feature('numCores')}, ...
        {'nMaxThreads',                 64}, ...
        ...                                                             % Plotting
        {'savepath',                    []}, ...                        %   field for saving the output path
        {'nFinePoints',                 300}, ...                       %   number of fine time points for plotting
        {'par_close_to_bound',          0.01}, ...                      %   notify if parameter within 1% of bound (relative to ub-lb)
        {'nfine_dr_method',             'pchip'}, ...                   %   spline
        {'nfine_dr_plot',               0}, ...                         %   use value > 10 to enable smoothing of DR curves
        {'plot_x_collected',            false}, ...                     %   0 = show seperate subplot for species and inputs, 1 = all in one
        {'ploterrors',                  0}, ...                         %   plotting options of error bars: 0=like fitted (error bar if yExpStd available, error band otherwise), 1=only error bars,  2=only error model as error band,  -1=confidence bands, -3: no errors
        {'showFitting',                 0}, ...                         %   Show the fitting process in real time
        {'showLegends',                 true}, ...                      %   Show legends in plots
        {'useSuptitle',                 false}, ...
        ...                                                             % Stochastic simulation
        {'ssa_min_tau',                 1e-3}, ...                      
        {'ssa_runs',                    1}, ...
        ...                                                             % Fit error handling
        {'fiterrors',                   0}, ...                         %   Fit error models?
        {'fiterrors_correction',        1}, ...                         %   Field for storing the Bessel-like error correction
        {'fiterrors_correction_warning',false}, ...                     %   Field for storing whether the user has been warned of the disabled Bessel-like error correction
        {'useFitErrorCorrection',       true}, ...                      %   Use Bessel-like correction when fitting error parameters
        {'useFitErrorMatrix',           false}, ...
        {'add_c',                       50},...                         % additive constant required in arCalcRes.m for lsqnonlin in case of error-fitting
        ...                                                             % Sampling
        {'useLHS',                      false}, ...                     %   When sampling random parameters use Latin Hypercube Sampling    
        ...                                                             % Optimization options
        {'useSensis',                   true}, ...                      %   Use sensitivities
        {'sensiSkip',                   false}, ...                     %   Skip sensitivities during fitting when only func is requested (speed-up for some optimizers)
        {'useJacobian',                 true}, ...                      %   Use Jacobian
        {'useSparseJac',                false}, ...                     %   Use Sparse Jacobian
        {'useSensiRHS',                 true}, ...                      %   Use sensitivities of RHS during simulation
        {'atolV',                       false}, ...                     %   Observation scaled tolerances
        {'atolV_Sens',                  false}, ...                     %   Sensi tolerances?
        {'optimizer',                   1}, ...                         %   Default optimizer
        {'optimizers',                  {'lsqnonlin', 'fmincon', 'PSO', 'STRSCNE', 'arNLS', 'fmincon_as_lsq', 'arNLS_SR1',...
                                         'NL2SOL','TRESNEI','Ceres', 'lsqnonlin_repeated', 'fminsearchbnd', 'patternsearch',...
                                         'patternsearch_hybrid', 'particleswarm', 'simulannealbnd', 'geneticalgorithm', 'lsqnonlinHeuristics'} }, ...
        ...                                                             % CVODES settings
        {'atol',                        1e-6}, ...                      %   Absolute tolerance
        {'rtol',                        1e-6}, ...                      %   Relative tolerance
        {'maxtol',                      1e-10}, ...                     %   Minimal tolerance used in regulating algorithms
        {'useTolTrustPar',              0}, ...                         %   Use parallel shrinkage of integration tolerances with optimizer trust region
        {'useTolSwitching',             0}, ...                         %   Use switch to strict integrator tolerances if trust region step is rejected
        {'maxsteps',                    1000}, ...                      %   Maximum number of steps before timeout
        {'maxstepsize',                 1e6}, ...                       %   Maximum stepsize
        {'useEvents',                   0}, ...                         %   Use event system
        {'useMS',                       0}, ...                         %   Use multiple shooting (DEPRECATED)
        {'nCVRestart',                  NaN}, ...                        %   Maximum number of automatic restarts
        ...                                                             % Simulation based equilibration settings
        {'init_eq_step',                100.0}, ...                     %   Simulation time of initial equilibration attempt
        {'eq_tol',                      1e-8}, ...                      %   Value below which all components of dxdt have to fall to be considered equilibrated
        {'eq_rtol',                     1e-8}, ...                      %   Relative value below which all components of dxdt have to fall to be considered equilibrated
        {'max_eq_steps',                20}, ...                        %   Maximum number of times the equilibration time is extended
        {'eq_step_factor',              5}, ...                         %   Factor by which the equilibration time is extended when dxdt isn't below eq_tol
        ...                                                             % Rootfinding based equilibration settings
        {'rootfinding',                 0},...                          %   Determine steady states by rootfinding rather than simulation
        ...                                                             % Constraint based steady states
        {'steady_state_constraint',     1}, ...                         %   Enable system
        ...
        {'instantaneous_termination',   1}, ...                    	% Poll utIsInterruptPending() to respond to CTRL+C
        {'no_optimization',             0}, ...                         % Disable compiler optimization
        };
      
    % Apply the default general settings where no fields are present
    ar.config = validateFields(ar.config, defaults, 'config');
    
    % PPL options
    pplDefaults = { ...
        {'alpha_level',           0.05}, ... 
        {'qLog10',                   0}, ...
        };
            
    % Apply the default PPL settings where no fields are present
    ar.ppl = validateFields(ar.ppl, pplDefaults, 'ppl');
    
    % PLE options
    pleDefaults = { ...
        {'alpha',                 0.05}, ... 
        {'ndof',                     1}, ...
        };
            
    % Apply the default PPL settings where no fields are present
    ar.ple = validateFields(ar.ple, pleDefaults, 'ple');
    
    if ~isfield( ar.config, 'useNewPlots' )
        if (isunix)
            ar.config.useNewPlots = false;
        else
            ar.config.useNewPlots = true;
        end
    end
    
    if ~isfield( ar.config, 'optim' )
        ar.config.optim = optimset('lsqnonlin');
        ar.config.optim.Display = 'off';
        ar.config.optim.TolFun = 0;
        ar.config.optim.TolX = 1e-6;
        ar.config.optim.MaxIter = 1000;
    end
    
    % Ceres optimization options
    if ~isfield( ar.config, 'optimceres' )
        ar.config.optimceres.TrustRegionStrategyType = 1;
        ar.config.optimceres.TrustRegionStrategyTypes = {'DOGLEG', 'LEVENBERG_MARQUARDT'};
        ar.config.optimceres.DoglegType = 1;
        ar.config.optimceres.DoglegTypes = {'SUBSPACE_DOGLEG', 'TRADITIONAL_DOGLEG'};
        ar.config.optimceres.LossFunctionType = 1;
        ar.config.optimceres.LossFunctions = {'None', 'TrivialLoss', 'HuberLoss', 'SoftLOneLoss', 'CauchyLoss', 'ArctanLoss'};
        ar.config.optimceres.LossFunctionVar = 1;
        ar.config.optimceres.TolFun = 0;
        ar.config.optimceres.TolX = 1e-6;
        ar.config.optimceres.TolGradient = 0;
        ar.config.optimceres.MaxIter = 1000;
        ar.config.optimceres.useNonmonotonicSteps = false;
        ar.config.optimceres.maxConsecutiveNonmonotonicSteps = 5;       
        ar.config.optimceres.maxSolverTimeInSeconds = 1e6;        
        ar.config.optimceres.NumThreads = 1;
        ar.config.optimceres.NumLinearSolverThreads = 1;       
        ar.config.optimceres.InitialTrustRegionRadius = 1e4;
        ar.config.optimceres.MaxTrustRegionRadius = 1e16;
        ar.config.optimceres.MinTrustRegionRadius = 1e-32;
        ar.config.optimceres.MinRelativeDecrease = 1e-4;
        ar.config.optimceres.MinLMDiagonal = 1e6;
        ar.config.optimceres.MaxLMDiagonal = 1e32;
        ar.config.optimceres.MaxNumConsecutiveInvalidSteps = 10;
        ar.config.optimceres.JacobiScaling = true;
        ar.config.optimceres.useInnerIterations = false;
        ar.config.optimceres.InnerIterationTolerance = 1e-3;
        ar.config.optimceres.LinearSolverType = 1;
        ar.config.optimceres.LinearSolvers = {'DENSE_QR','DENSE_NORMAL_CHOLESKY', 'CGNR', 'DENSE_SCHUR', 'SPARSE_SCHUR', 'ITERATIVE_SCHUR'};
        ar.config.optimceres.printLevel = 0;
	end
    
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
    ar.info.arsimucalc_flags{15} = 'CVodeSetMaxNumSteps()';
    ar.info.arsimucalc_flags{16} = 'CVodeReInit()';
    ar.info.arsimucalc_flags{17} = 'CVodeSensReInit()';
    ar.info.arsimucalc_flags{18} = 'malloc EventData';
    ar.info.arsimucalc_flags{19} = 'CVodeSetMaxNumSteps()';
    ar.info.arsimucalc_flags{20} = sprintf('equilibration. Failed to meet tolerance.\nDoes the system have a steady state? Failure occurred');
    ar.info.arsimucalc_flags{21} = sprintf('initial condition override.\nInitial condition override vector has the wrong size');

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

    ar.info.arFormatVersion  = arFormatVersion;
    
function str = validateFields(str, fields, fieldname)

    for jf = 1 : size( fields, 2 )
        if ( ~isfield( str, fields{jf}{1} ) )
            str.(fields{jf}{1}) = fields{jf}{2};
            arFprintf( 3, '  Added %s.%s\n', fieldname, fields{jf}{1} );
        end
    end
    
function str = checkForField(str, field)
    if ( ~isfield( str, field ) )
        str.(field) = struct;
        arFprintf( 3, 'Added substructure %s\n', field );
    end
