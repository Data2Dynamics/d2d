
function RealisticDesign_D2D(plog, p2)

% Create a realistic design with Observables, Time Points, Simulated Data and 10% Errors
% Biomodel and transient model have to be created/compiled with Setup.m first
% ar = arStruct of model from biomodels database

% RealisticDesign_D2D:
%1) Define a realistic sample of Observables out of given model states (direct, scaled, compound measurements possible)
%2) Fit a Transient Function to the observables dynamic (TransientFct.pdf)
%3) Calculate realistic Time Points with time parameters of Transient Function (SuggestedTimePoints.pdf)
%4) Simulate data with mean from model dyanmics and std deviation of 10% (Data.xls)
%5) Save Realistic Design in ar Struct (workspace.mat)
%
% Inputs:
% ar = arStruct of model from biomodels database
% plog = Do you want to log the parameters? Then set true.(Default: false)
%        But watch out! Logging parameters, e.g. if ar.p = 0 can cause
%        differences to the original model
% changeT = If you want to change time scale --> Unrealistic model

global ar

% Check existence of variables
if ~exist('plog','var') || isempty('plog')
    fprintf('No information about logging given. No log applied.\n')
    plog = false;
end
if plog == true
    warning('Watch out! Logarithm of parameters can cause differences to original model.');
end
if ~exist('p2','var') || isempty('p2')
    p2 = false;
end
if ~exist('changeT','var')
    fprintf('Realistic Design will be created.\n')
    changeT = [];
else
    fprintf(['Unrealistic Design with ' changeT '*(time interval) will be created.\n'])
end

% remove files from previous run, create info folder

if exist([pwd '/RealisticDesign'],'dir')
	rmdir([pwd '/RealisticDesign'],'s');
end
if ~exist([pwd '/RealisticDesign'],'dir')
    mkdir('RealisticDesign');
end


%% Load biomodel
if ~exist('ar','var') || isempty(ar)
    if exist('Results/Biomodel','dir')
        arLoad('Biomodel')
    else
        fprintf('First load biomodel')
        if exist('Results','dir')
            arLoad
        else
            error 'No Model found. Did you compile one?\n'
        end
    end
end


%% State Observables
Observables;  %% Dice Observables from States

%% Convert order of magnitude
%Magnitude_Conversion;

%% Fit Transient
arModel = arDeepCopy(ar);  % Need ar of Transient, so remember Biomodel
taus = FitTransient;
ar = arDeepCopy(arModel);

%% Get Time Points
TimePoints_adddt(taus);

%% Write data def with new Obs def (+offset/scale)
WriteDataDef(plog,p2);
%arSave('Real')
%% Simulate Data Points
RealisticSimu(plog,p2);

