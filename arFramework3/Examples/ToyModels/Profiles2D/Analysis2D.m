%% Calculate Parameter Profiles

arPLEInit([],[],[],false);
ar.ple.relchi2stepincrease = 0.05 * ones(size(ar.ple.p));
    % More accurate Profile Sampling
ple([],100);
pleSmooth;
ar.config.savepath = '.\Results\Profiles';
ar.config.backup_modelAndData = false;
arSave('current',[],false);

%% Calculate Validation Profile for the experimental condition of interest

InitVPL(1,1,2,40,0.2);
    % Validation profile for observable = 2 (B), t = 40 and sigma = 0.2
ar.vpl.config.chi2dif_max = 0.2;
    % More accurate profile sampling
ar.vpl.config.maxsteps = 200;
VPL;
vplSmooth;

% vplPlot

% Output is saved in ar.vpl

%% Calculate 2D-Profiles

%%% Generating the 2D-Profile

Init2DPL(10); 
    % Number of 2d profiles in each direction from optimum
ar.ple2d.config.gen2d.weakness = 10; 
    % Weak quadratic prior if validation profile is non-identifiable
gen2d(2);
    % Generate the 2D-Profile for parameter 2
autofix2dprofile;
    % Automatically diagnose and fix discontinuities
% Output is saved in ar.ple2d   
   
%%% Analyzing the 2D-Profile

bounds2d; 
    % Calculates 95% confidence bounds for parameters at different
    % predictions (also saves it in ar.ple2d for plotting)
s = score2d;
    % Calculates the expected average parameter profile width
    
% Plot2DPL(2);