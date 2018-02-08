% Often, the dynamics of the model should be in a steady state at the beginning 
% of a simulation. This reflects the biological situation of cells being in a 
% resting (i.e. unstimulated) state. D2D offers three approaches to deal 
% with steady states; each of which have their own advantages and drawbacks
%
% ANALYTICAL STEADY STATES
%   Analytical steady states via parameter transformations in the CONDITIONS field.
%     Advantages: 
%       - No simulation required. Reduction in the number of parameters.
%     Disadvantages: 
%       - Only tractable for relatively simple steady states. 
%         Up to the user to correctly determine the steady state equations
%     Method:
%       Analytical steady states can be imposed by simply enforcing them in
%       the CONDITIONS section of a data def file.
%
% SIMULATED STEADY STATES
%   Steady state imposition by equilibrating the system to steady state.
%     Advantages: 
%       - No additional term in optimization.
%       - Reduction in the number of parameters.
%     Disadvantages: 
%       - Pre-simulation required. 
%       - Integration failure when steady state does not exist.
%     Method:
%       Pre-equilibration can be applied by using arSteadyState.
%       Invoke: help arSteadyState for more information on how to set this up
%
% CONSTRAINT BASED STEADY STATES
%   Steady state imposition by means of penalized optimization.
%     Advantages: 
%       - Optimized initial condition parameters correspond to a (near) 
%         steady state and it is also possible to consider systems 'near' steady state.
%     Disadvantages:
%       - Constraint is included in profiling and optimization, which can 
%         trade off penalty and fidelity to the data. 
%       - Have to critically evaluate what the parameter profiles really mean 
%         when introducing such a constraint. Since constraints are ratio-metric, 
%         states near zero are penalized much more strongly in absolute terms.
%       - No clear statistical method for choosing the penalty weight.
%     Method:
%       For a specific model jm and condition jc set:
%         ar.model(jm).condition(jc).qSteadyState(:) = 1;
%       Subsequently set the inverse weight:
%         ar.model(jm).condition(jc).stdSteadyState(:) = std_val;
%       Note: You can find the condition using arShowDataConditionStructure.
%

help arHelpSteadyState