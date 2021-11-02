function Init2DPL(npred_chi,del)
% Init2DPL(npred_chi,del)
%
% Initializes the struct for the 2D-Profile Likelihood. Initialization of
% the method should be from the global optimum. Entries from a previous
% 2d-calculation are removed when calling this function.
%
% npred_chi: Number of prediction steps in each direction of the minimum
%            for which parameters profiles are calculated. [40]
% del:       This is the stepsize in alpha-space as an alternative to
%            specifying npred_chi. This will be overwritten if npred_chi is
%            specified.  
%
% See also: gen2d, ple, VPL

global ar

%% Initial checks:

% Clear entries from a previous 2d calculation:
if isfield(ar,'ple2d')
    ar = rmfield(ar,'ple2d');
end

% Match parameter profile confidence level with the level used for the
% original calculation
try
    ar.ple2d.config.bounds.alpha_bounds = 1-ar.ple.alpha_level;
catch
    disp(['WARNING Init2DPL: Calculate respective Parameter Profile at',...
        ' a specified confidence level first'])
    return
end

if isfield(ar,'vpl')
    alpha_vpl = 1-ar.vpl.config.alpha; 
    max_vpl_level = 1.1*icdf('chi2',alpha_vpl,1);
    % Specifies the maximal VPL-value up to which parameter profiles should 
    % be calculated.
else
    disp('WARNING Init2DPL: Calculate the validation profile first by calling VPL')
    return
end

%% Necessary settings for gen2d:
% These setting implicitly determine at which predictions parameter
% profiles are calculated by constructing a set of corresponding validation 
% profile levels. Steps are constant steps in alpha-space , i.e. constant 
% changes of confidence level if validation profile is sufficiently asymptotic. 

% Find number of steps and step sizes in alpha-space:
if(~exist('npred_chi','var') || isempty(npred_chi))
    if(~exist('del','var') || isempty(del))
        npred_chi = 40;
        del = cdf('chi2',max_vpl_level,1)/npred_chi;
    else
        npred_chi = ceil(cdf('chi2',max_vpl_level,1)/del);
    end
else
    del = cdf('chi2',max_vpl_level,1)/npred_chi;
end
ar.ple2d.config.gen2d.npred_chi = npred_chi;
ar.ple2d.config.gen2d.del = del;
ar.ple2d.config.gen2d.levels = icdf('chi2',del*(1:npred_chi),1);
% These are the levels at which parameter profiles are calculated

ar.ple2d.config.gen2d.weakness = 15;
% weakness times sigma is the width of the quadratic prior

% Options for ple which is called in gen2d:
ar.ple2d.config.gen2d.plmode = 2;
ar.ple2d.config.gen2d.npar = 100;
ar.ple2d.config.gen2d.relchi = 0.1;

%% Options for other 2D-functions:

ar.ple2d.config.bounds.save_mode = 2;
% if save_mode is 1 bounds will not be saved into ar struct when calling bounds2d,
% it is 2 otherwise. Always leave it on 2, option 1 is only used in score2d and
% it is enabled automatically.
ar.ple2d.config.bounds.vpl_mode = 1;
% vpl_mode determines whether grid interpolated values (vpl_mode = 2) or actual
% parameter and validation profile values (vpl_mode = 1) should be used.
ar.ple2d.config.bounds.nvpl = 20;
% nvpl is a guiding value for the number of points which are used for
% averaging in each validation direction
ar.ple2d.config.bounds.alphavpl_max = alpha_vpl;
ar.ple2d.config.bounds.alphavpl_min = 0.1;
% alphamax/alphamin are guide values between which alpha values prediction
% values should be regarded. This does not concern the profile sampling in
% gen2d.

% Plot options:
ar.ple2d.config.plot.thicknessscatter = 60;
ar.ple2d.config.plot.ncontourlines = 40;
ar.ple2d.config.plot.xnodes = 800; % also used for smooth2d

% Autofix options:
ar.ple2d.config.autofix.excess_factor = 1.2;
    % Confidence Threshold Modifier (SmoothingPoints2d_Excess, extend2dprofile)
ar.ple2d.config.autofix.eps_inside = 0.05; 
    % Threshold for inside profile jumps (SmoothingPoints2d_Inside)
ar.ple2d.config.autofix.eps_outside = 0.05;
    % Threshold for outside profile jumps (SmoothingPoints2d_Bounds)
ar.ple2d.config.autofix.eps_sample = 0.05;
    % Threshold for inside sampling errors (rmInsideJumps2d)
ar.ple2d.config.autofix.niter = 3;
    % Number of fixing iterations before current method is skipped
ar.ple2d.config.autofix.jumpfactor = 2;
    % Lower values will detect more possible jumps
    % (SmoothingPoints2d_Excess, SmoothingPoints2d_Inside)
ar.ple2d.config.autofix.scalefactor = 10^-4;
    % Determines minimum size of jumps which can be detected (multiplied by
    % range of parameter profile) (SmoothingPoints2d_Excess);

end






