function model = StrucIDreadmodel(modelname)
% modelname refers to a text file that contains the model definition in a
% specific format (see example file)
% The result Model is a structure with fields

% run ADiMAT startup script for automatic differentiation. This is to make
% sure that the ADiMAT routines are available when running the sensitivity
% computations
if (~exist('admDiffComplex.m','file') && ~isdeployed)
    %run('./adimat/ADiMat_startup');
    ADiMat_startup;
end

%% 1. ENTER MODEL _ GO AND FETCH FILE TO ATTACH

if ~ischar(modelname)
    error('Error. Argument to this function must be a string of characters.');
else
    filename = append(modelname,'.txt');
    delimiterIn = '!'; %delimiter used in code to detect different sections
    tree = importdata(filename,delimiterIn);
end

close all; % close all figures if these are still open
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[algebraicEqns,ODEs,stateNames,stateValues,ICStateEqns,parNames,parValues,UnknownVarIndex,...
    inputNames,inputEqns,outputNames,outputEqns,~] = process_model_file(tree);

% app.resultsFolder=pathname;

% vAlgebraText=algebraicEqns;
% XdotText=ODEs;
% YoutText=outputEqns;
% UinText=inputEqns;
% app.kReducedEditField.Value=numel(ODEs);

% Build functions of the algebraic equations, ODEs, Yout, Inputs, and x0.
% The function handles + definitions will be passed on as outputs/fields to
% the model structure variable with name 'model'

% algebraic relations
vAlgebra = str2func(['@(x,p,u)[',strjoin(algebraicEqns),']']);

% ODEs or dynamics f(t,x,u,p,v)
Xdot = str2func(['@(t,x,u,p,v)[',strjoin(ODEs),']']);

% Initial conditions
x0State = str2func(['@(p)[',strjoin(ICStateEqns),']']);

% output equations (or sensors available)
Yout = str2func(['@(t,x,u,p,v)[',strjoin(outputEqns),']']);

% inputs/stimuli for the dynamics
Uin = str2func(['@(t,p)[',strjoin(inputEqns),']']);

% set default variables for the identifiability analysis
integrator = @ode45; % default integration method
Tf=3; % final time integration
JacobiComplex = true; % default is to use complex variables for the integration of the sensitivity eqns
nTraject = 1; % integrate only one trajectory by default
FixedPar = ~isnan(parValues); np=numel(parValues);
FixedState = ~isnan(stateValues); nx=numel(stateValues);
FixedVar = cat(1,FixedPar,FixedState);

if any(FixedVar) % check if any elements are non-zero
    index=1:(np+nx); FixedVarIndex=index(FixedVar);
else
    FixedVarIndex=[];
end

model = struct('fDyn',Xdot,'U',Uin,'vAlg',vAlgebra,'x0',x0State,'Y',Yout,...
    'stateNames',{stateNames},'stateValues',stateValues,'parNames',{parNames},'parValues',parValues,...
    'inputNames',{inputNames},'outputNames',{outputNames},'UnknownVarIndex',UnknownVarIndex,'integrator',integrator,...
    'Tf',Tf,'JacobiComplex',JacobiComplex,'nTraject',nTraject,'FixedVarIndex',FixedVarIndex);

