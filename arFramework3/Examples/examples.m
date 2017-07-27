% In this folder you will find several examples illustrating the use of
% various parts of the D2D framework. It is not recommended to edit these
% examples, but work with copies of them if you want to base future work on
% examples from this directory, since some of them are used for testing the
% framework (see doTests.m).
% 
% The root contains finished models which resulted in various papers:
% Bachmann_MSB2011, Becker_Science2010, Beer_MolBiosyst2014,
% Boehm_JProteomeRes2014, Bruno_Carotines_JExpBio2016,
% Merkle_JAK2STAT5_PCB2016, Oxygen_Metabolic_Rates2016,
% Raia_CancerResearch2011, Schwen_InsulinMouseHepatocytes_PlosOne2014,
% Swameye_PNAS2003, Toensing_InfectiousDisease_BoardingSchool_2017,
% Toensing_InfectiousDisease_Zika2017
% as well as the models used in the dream challenge: Dream6, Dream6_L1
%
% The instructive examples treat various advanced topics in D2D and are 
% structured as follows:
%
% ToyModels/
%   .SubSensitivities
%       This toymodel demonstrates how to use subsensitivities. For large
%       models, evaluating the sensitivity equations can become
%       prohibitively expensive. In such cases, one might opt for
%       optimizing only a subset of the parameters (for example the Vmaxes
%       in a metabolic model). In this case, it is wise to activate the
%       subsensitivity system. See this example for an example how to
%       activate this system.
%   Input_Tests
%       This toymodel demonstrates various input functions (bolus
%       injection, step function, double step function, smooth step
%       functions)
%   Step_Estimation
%       This toymodel demonstrates how to estimate parameters in input
%       functions which behave in a (near) stepwise manner. It also shows
%       how to use smooth step functions which have a varying smoothness
%       parameter (and are therefore better at estimation the location of
%       the step)
%   Advanced_Events
%       This folder contains various examples demonstrating the event
%       system.
%       - EquilibrationExample shows how to numerically equilibrate multiple 
%         conditions and use them as initial values for subsequent simulations
%       - eventExample shows how to add run-time events which modify a state
%         instantaneously at run-time (without recompilation)
%       - conservedExample, shows how to use the rootfinding procedures
%         which avoid numerically simulating the steady state (it depends on
%         your model whether this is slower or faster than simulated steady
%         states)
%   Intercondition_constraints
%       -  This example shows how to use regularization between conditions.
%          It shows how to add a constraint between two conditions such
%          that two conditions will get the same simulation if they do not
%          require to be different for the data.
%   LongSplines
%       - This example shows how to string together splines in order to
%         use a long driving input with more than 10 points.
%   TurboSplines
%       - This example shows how to activate the turbosplines system which
%         caches the spline coefficients for use in the RHS. This system
%         is much faster to evaluate but (slightly) slower to compile.
%   Volume_Estimation
%       - This example shows an example how to estimate the volume of
%         compartments along with the rest of the problem
%   Splines
%       - This example shows the use of various types of cubic splines.
%   SteadyStateModel
%       - This example shows how to fit a steady state model and skip the
%         dynamic simulation.
%   Flux_Estimation
%       - This example shows how to use fluxes in the obsevation function
%   Bolus_Injection_Test
%       - This model shows how to simulate a bolus injection
%   Washing_and_Injection_Test
%       - This model shows how to simulate a bolus injection and washing
%         experiment

help examples