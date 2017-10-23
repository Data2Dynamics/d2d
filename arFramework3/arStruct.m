%
% The ar struct contains most of the data used for all model analyses
%
%   Model - Model
%     - res                 Residual
%     - sres                Sensitivity of the residual (here the observables and parameters are specified in whatever space they are fitted (log or not))
%     - constr              Constraint based residuals
%     - chi2                Sum of squared residuals (when fitting errors, this is *not* the objective function)
%                             Given by sum( ar.model.data.res^2 )
%     - chi2err             Sum of error residuals
%                             Given by sum(ar.model.data.reserr.^2) + ar.config.add_c * sum(ar.model.data.reserr>0)
%                             These come from the normalization constant of the likelihood which needs to be taken into account when estimating sigma.
%                             Note that the add_c term comes from the fact that we need to provide a residual to lsqnonlin. Since the normalization term is proportional to
%                             sqrt(ln(sigma)+add_c), the add_c is used to ensure positivy under the sqrt. It is later subtracted again when computing the objective function.
%     - chi2constr          Sum of squared constraint residuals.
%                             Given by sum(ar.constr.^2)
%     - chi2fit             Sum of squared residuals and error residuals (this is typically the objective function when there are no constraints).
%                             Given by ar.chi2 + ar.chi2err
%                             Hint: the log likelihood can be computed as 2*ar.ndata*log(sqrt(2*pi)) + ar.chi2fit;
%     - (chi2s)             Present after running an arFits or LHS. List of objective functions with NaN's for ones that did not complete.
%     - (ps)                Present after running an arFits or LHS. List of parameter vectors with NaN's for ones that did not complete
%     - pLabel              Parameter names
%     - qFit                Vector of logicals whether parameter is fixed (0) or fitted (1)
%     - qLog10              Vector of logicals whether parameter should be considered in logspace
%     - qDynamic            Vector of logicals which indicate whether parameter affects dynamics
%     - qError              Vector of logicals which indicates whether it is an error model parameter
%     - qInitial            Vector of logicals which indicates whether a parameter is an initial condition
%     - lb/ub               Lower and upper bounds for the parameters (note that depending on qLog10 they are in logspace or not)
%     - p                   Parameter values (note that depending on qLog10 they are in logspace or not)
%     - type                Type of prior (0=box, 1=normal, 2=uniform with normal bounds, 3=L1 prior)
%
%     - condition           Experimental conditions
%         dLink             Index to which data this condition belongs
%         tExp              Time points used in fitting
%         tFine             Time points for the fine simulation
%         tEvents           Event time points (lead to solver
%                           re-initialization)
%         modx_A-modsx_B    Part of the event system. These override the
%                           states and their sensitivities Ax+B.
%         xExpSimu          State trajectories for this condition
%         xFineSimu         State trajectories for this condition
%         uExpSimu          Input trajectories for this condition
%         uFineSimu         Input trajectories for this condition
%         vExpSimu          Flux trajectories for this condition
%         vFineSimu         Flux trajectories for this condition
%         zExpSimu          Derived variable trajectories for this condition
%         zFineSimu         Derived variable trajectories for this condition
%         pLink             Which external parameters are used in this
%                           condition?
%         sxExpSimu         Model state sensitivities for fitting [nExp,nState,nPars]
%                               These sensitivities are always specified in
%                               linear parameter space.
%         sxFineSimu        Fine model state sensitivities  [nFine,nStates,nPars]
%                               When parameters are being fitted in
%                               logspace, these sentivities are specified
%                               with respect to logspace.
%                               Note: dx/dlog10(p) = dx/dp * ln(10) * 10^log10(p)
%
%     - ss_condition        Steady state conditions used for equilibration
%                           (mostly analogous to condition)
%         src               Which condition to use as reference for the
%                           steady state equilibration.
%         ssStates          Which states to equilibrate and carry over.
%         ssLink            Which conditions to transfer the simulated
%                           steady state to.
%
%     - data                Data structure
%         y                 Names of the observables
%         cLink             Index to which conditions correspond to this data
%         tExp              Time points used in fitting
%         tFine             Time points for the fine simulation
%         yExpSimu          Observable trajectories for this data
%         yFineSimu         Observable fine grained trajectories
%         res               Residual vector corresponding to the data
%         sres              Sensitivities of the residual vector
%         reserr            Part of the residual vector that corresponds to
%                           the normalization constant of the likelihood
%                           sqrt(ln sigma + add_c). Note that a constant is added
%                           to ansure positivity.
%         sreserr           Sensitivities of reserr.
%         syExpSimu         Observable state sensitivities for fitting [nExp,nState,nPars]
%                               When fitting in logspace, these are specified in log10(y).
%                               Note: dlog10(y)/dp = (dy/dp) / (y*ln(10)) [This is performed in arSimuCalc.c]
%                               These sensitivities are specified in *linear parameter space*.
%         syFineSimu        Fine model state sensitivities  [nFine,nStates,nPars]
%                               When fitting in logspace, these are specified in log10(y).
%                               When parameters are fitted in logspace, these sentivities are specified in *logspace*.
%
%     - Plots               Contains information on how the plots are
%                           linked to the data
%
