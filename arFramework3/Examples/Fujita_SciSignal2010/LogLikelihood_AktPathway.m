function varargout = LogLikelihood_AktPathway(theta, D)
% Objective function for the Akt Pathway example
% 
%
% logLikelihood_AktPathway.m provides the log-likelihood and its gradient 
% for the model defined in AktPathway_syms.m
% 
% USAGE:
% [llh] = LogLikelihood_AktPathway(theta, amiData)
% [llh, sllh] = LogLikelihood_AktPathway(theta, amiData)
%
% Parameters:
%  theta: Model parameters 
%  amiData: an amidata object for the AMICI solver
%
% Return values:
%   varargout:
%     llh: Log-Likelihood, only the LogLikelihood will be returned, no 
%         sensitivity analysis is performed
%     sllh: Gradient of llh, The LogLikelihood and its gradient will be 
%         returned, first order adjoint sensitivity analysis is performed



%% Model Definition
% The ODE model is set up using the AMICI toolbox. To access the AMICI
% model setup, see AktPathway_syms.m



%% AMICI
% Setting the options for the AMICI solver
optionsAmici = amioption;
optionsAmici.atol = 1e-8;
optionsAmici.maxsteps = 2e5;
optionsAmici.interpType = 2;

logL = 0;
var = zeros(16);

if(nargout > 1)
    for i=1:length(D)
        optionsAmici.sensi = 1;
        optionsAmici.sensi_meth = 'adjoint';
        sol = simulate_AktPathway_model(D(i).t, theta, D(i).condition, D(i), optionsAmici);
        if(sol.status < 0)
            error('integration error');
        else
            logL = logL + sol.llh;
            var = var + sol.sllh ;
        end
    end
    varargout{1}=logL;
    varargout{2}=var;
else
    for i=1:length(D)
        optionsAmici.sensi = 0;
        sol = simulate_AktPathway_model(D(i).t, theta, D(i).condition, D(i), optionsAmici);
        if(sol.status<0)
            error('integration error');
        else
            logL = logL + sol.llh;
        end
    end
    varargout{1}=logL;
end

end

