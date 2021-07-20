function [method,obj_fun] = ssm_mix(i,nrun)

%--------------------------------------------------------------------------
% Define here the fraction of searches to perform with each method, and the
% corresponding objective function:
DHC       = 0.3; % Value between 0 and 1
DHC_f     = 'like_MEIGO_B2'; % Objective function to use with DHC
NL2SOL    = 0.2; % Value between 0 and 1
NL2SOL_f  = 'likelihood_MEIGO_B2_NL2SOL'; % Objective function to use with NL2SOL
FMINCON   = 0.5; % Value between 0 and 1
FMINCON_f = 'like_MEIGO_B2'; % Objective function to use with FMINCON 

%--------------------------------------------------------------------------
if (DHC+NL2SOL+FMINCON ~= 1.0)
    disp('Warning: the sum of the selected percentages of local methods is not equal to 100%. Default percentages will be assigned.')
    % Default values:
    DHC     = 0.3; 
    NL2SOL  = 0.2; 
    FMINCON = 0.5; 
end
if i <= round(nrun*DHC)
    method = 'dhc'
    obj_fun = DHC_f;
else if i <= ( round(nrun*DHC) + round(nrun*NL2SOL) )
        method = 'nl2sol'
        obj_fun = NL2SOL_f;
    else method = 'fmincon'
        obj_fun = FMINCON_f;
    end
end

end

