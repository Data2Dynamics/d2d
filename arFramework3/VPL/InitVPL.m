function InitVPL(sigma,chi2dif_max,stepmode,maxsteps)
% InitVPL(sigma,chi2dif_max,mode,maxsteps)
%
% Initializes the struct ar.vpl for the validation profile calculation.
%
% sigma:       Standard Error of validation data point        [mandatory user input]
% chi2dif_max: Guiding value to control chi2 stepsize         [0.2]
% stepmode:        Defines which step choice algorithm is employed [2]
%       -1: Basic step choice based on chi2 difference of last step
%       -2: Adaptive step choice based on a theoretical upper chi2 limit
% maxsteps:    Maximal number of steps in each direction      [100] 
%
% Separating this function from the main function VPL allows easy
% manipulation of algorithm options by changing values in the new struct.

global ar

%Remove previous validation profile calculation if it exists
if isfield(ar,'vpl')
    ar = rmfield(ar,'vpl');
end
if ~exist('sigma','var')
    disp('ERROR InitVPL: Specify standard deviation sigma.')
    return
end

%Default values if not specified
if (~exist('chi2dif_max','var') || isempty(chi2dif_max))
    chi2dif_max = 0.2;
end
if (~exist('stepmode','var') || isempty(mode))
    stepmode = 2;
end
if (~exist('maxsteps','var') || isempty(maxsteps))
    maxsteps = 100;
end

%Technical algorithm inputs and magic factors
ar.vpl.general.sigma = sigma; 
%Sigma is defined in this function so that range and stepsizes can
%already be set accordingly
ar.vpl.config.chi2dif_max = chi2dif_max; 
ar.vpl.config.chi2dif_min = chi2dif_max/4; %Not actual minimimal value
ar.vpl.config.chi2max = 1.2*icdf('chi2',0.95,1); %Maximal chi2 profile value
ar.vpl.config.maxsteps = maxsteps; 
ar.vpl.config.stepfactor = 1.5; %Step size adaption factor
ar.vpl.config.maxrange = 50*sigma; %Maximal absolute range of profile in each direction
ar.vpl.config.maxstepsize = 5*sigma; %Maximal absolute value of step
ar.vpl.config.minstepsize = 5*sigma/1000; %Minimal absolute value of step
ar.vpl.config.firststep = sigma/20; %first step from the optimum 
%Optional flags:
ar.vpl.config.sensi = true; %integrate with sensitivities after fitting
ar.vpl.config.showCalculation = true; %print intermediate calculation results 
ar.vpl.config.prediction = false; %true: validation profile range is set by
                                  %         the prediction profile values

%Step choice function
if stepmode == 1
    ar.vpl.config.step_fkt = @vplStepPrevious; 
elseif stepmode == 2
    ar.vpl.config.step_fkt = @vplStepDynamic;
end

end

