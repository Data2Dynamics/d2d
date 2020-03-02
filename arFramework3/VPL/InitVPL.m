function InitVPL(m,d,idpred,tpred,sigma,stepmode)
% InitVPL(m,d,idpred,tpred,sigma,stepmode)
%
% m:           Model index
% d:           Data index
% idpred:      Index of prediction of interest
% tpred:       Prediction time of interest
% sigma:       Standard deviation of data point [default: simulated value]
% stepmode:    Defines which step choice algorithm is employed [2]
%       -1: Basic step choice based on chi2 difference of last step
%       -2: Adaptive step choice based on a theoretical upper chi2 limit
%Setu
% Initializes the struct ar.vpl for the validation profile calculation.
% Separating this function from the main function VPL allows easy
% manipulation of algorithm options by changing values in the new struct.
%
% See also: VPL

global ar

%Remove previous validation profile calculation if it exists
if isfield(ar,'vpl')
    ar = rmfield(ar,'vpl');
end

try
    ar.vpl.general.m = m;
    ar.vpl.general.d = d;
    ar.vpl.general.idpred = idpred;
    ar.vpl.general.tpred = tpred;
catch
    disp('ERROR InitVPL: Mandatory input arguments ill-specified')
    return
end

if (~exist('stepmode','var')) || (isempty(stepmode))
    stepmode = 2;
end
% Simulate sigma if not specified:
if (~exist('sigma','var')) || (isempty(sigma))
    if ar.config.fiterrors == -1
        fprintf(['\n ERROR InitVPL: Only experimental errors are used',...
            '(ar.config.fiterrors = -1) \n','  Please specify sigma. \n']);
        return
    end
    disp('Standard deviation sigma not specified:. Simulating...')
    
    try
        arAddToData(m,d,idpred,tpred,0,2,1); %0 is a data point placeholder
        arLink(true); %Add Data and Link separately to suppress output
        arSimu(true,false); %This simulates the data errors
        sigma = ar.model(m).data(d).ystdExpSimu(...
            size(ar.model(m).data(d).yExp,1),idpred);
    catch exception
        fprintf(['\n ERROR InitVPL: Resetting ar struct. Sigma could not be simulated. \n',...
            'Error message: %s \n Line: %s \n'],...
            exception.message, sprintf('%i, ',exception.stack.line));
        arRemoveData(1,1,2,1); 
        arLink(true);
        return
    end
    
    fprintf('\n Simulated sigma = %0.4g \n',sigma);
    arRemoveData(1,1,2,1); %Removes data points added by arAddToData
    arLink(true);
end
ar.vpl.general.sigma = sigma;

%Technical algorithm inputs and magic factors
ar.vpl.config.chi2dif_max = 0.2; 
ar.vpl.config.chi2dif_min = 0.05; %Not actual minimimal value
ar.vpl.config.chi2max = 1.2*icdf('chi2',0.95,1); %Maximal chi2 profile value
ar.vpl.config.maxsteps = 100; 
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

