
function arRealisticDesign(realisticname,cond,plog, p2)

% Create a realistic design with Observables, Time Points, Simulated Data and 10% Errors
% Biomodel and transient model have to be created/compiled with Setup.m first
% ar = arStruct of model from biomodels database

% RealisticDesign_D2D:
%1) Define a realistic sample of Observables out of given model states (direct, scaled, compound measurements possible)
%2) Fit a Transient Function to the observables dynamic (TransientFit_*.png)
%3) Calculate realistic Time Points with time parameters of Transient Function (TimePoints.txt)
%4) Simulate data with mean from model dyanmics and std deviation (RealisticData.xls)
%5) Save Realistic Design in ar Struct (workspace.mat)
%
% Inputs:
% realisticname = string for naming the ar of realistic simulation   ['Realistic']
% cond = condition number                                       [1]
%           realistic design just works for one condition, 
%           biomodels usually come with just one condition,
%           if your model has multiple conditions, consider
%               idc = FindMostStatesInCondition;
%               idd = ar.model.condition(idc).dLink;
%               arQFit('data',0)
%               arQFit('data',1,idd)
% plog = Do you want to log the parameters?                     [true]
%        But watch out! Logging parameters, e.g. if ar.p = 0 can cause
%        differences to the original model
% p2 = Parameter bounds are set to [p-2,p+2] if true            [true]

global ar

% Check existence of variables
if ~exist('realisticname','var') || isempty(realisticname)
    realisticname = 'Realistic';
end
if ~exist('cond','var') || isempty('cond')
    cond=1;
end
if ~exist('plog','var') || isempty('plog')
    fprintf('No information about logging given. Logarithm of observables is drawn realisticly. \n')
    plog = true;
end
if plog == true
    warning('Watch out! Logarithm of intensities can cause differences to original model. If this is not what you want, set arRealisticDesign([],false). \n');
end
if ~exist('p2','var') || isempty('p2')
    p2 = true;
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
    if exist(['Results' filesep 'Biomodel'],'dir')
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

%% Set Observables
if ~isfield(ar.model,'data')
    arSetObservables(cond);  %% Dice Observables from States
    dataname = [];
else
    warning('%s\n','Observables are NOT drawn. Existing Observables rest.')
    dataname = ar.model.data.name;
end

%% Convert order of magnitude
convertt = arMagnitudeConversion;

%% Fit Transient
arModel = arDeepCopy(ar);  % Need ar for Transient, so remember Biomodel
taus = arTransientPars;
ar = arDeepCopy(arModel);

%% Time Points
arSetTimePoints(taus,length(ar.model.x),length(ar.model.data.y),convertt);

%% Compile & Simu Data
arRealisticSimu(ar.model.name,realisticname,plog,p2,dataname);
    