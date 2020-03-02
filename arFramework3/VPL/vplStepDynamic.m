function [delz] = vplStepDynamic(z_old,chi2dif,delz,sigma,gen_struct)
% delz = vplStepDynamic(z_old,chi2dif,delz,gen_struct)
%
% Chooses validation data step delz. Dynamic step size adaption by use of a
% theoretical upper limit as a test-value.
%
% z_old:      Most recent valdiation data point
% chi2dif:    Most recent actual chi2 difference (i.e. after fitting)
% delz:       Most recently chosen stepsize
% gen_struct: Various necessary variables from upper level function

%From the higher level function: 
pred_old          = gen_struct.temp.pred;
delii             = gen_struct.temp.delii;
%Magic factors:
chi2dif_max       = gen_struct.temp.chi2dif_max;
maxstepsize       = gen_struct.temp.maxstepsize;
minstepsize       = gen_struct.temp.minstepsize;
stepfactor        = gen_struct.temp.stepfactor;

%Theoretical upper limit for chi2 change:
chi2limit = (delii*delz/(sigma^2))*...
    (2*(z_old-pred_old)+delii*delz);

%Step size adaption chosen in a way that puts chi2limit as closely to the
%theoretical maximum as possible
if chi2limit > chi2dif_max
    while chi2limit > chi2dif_max
        if delz/stepfactor > minstepsize
            delz = delz/stepfactor;
        else
            delz = minstepsize;          
            break
        end
        chi2limit = (delii*delz/(sigma^2))*...
            (2*(z_old-pred_old)+delii*delz);
    end
else
    while chi2limit <= chi2dif_max
        if delz*stepfactor < maxstepsize
            delz = delz*stepfactor;
        else
            delz = maxstepsize;  
            delz = delz*stepfactor; %compensates operation directly after while loop
            break
        end
        chi2limit = (delii*delz/(sigma^2))*...
            (2*(z_old-pred_old)+delii*delz);
    end
    delz = delz/stepfactor; %Puts chi2limit back under the threshold
end

%Be more aggressive if last actual observed chi2-difference is too small.
if ~((delz == maxstepsize) || (delz == minstepsize)) && ...
        (abs(chi2dif) < gen_struct.temp.chi2dif_min) &&...
        (gen_struct.temp.correction_switch)
        delz = stepfactor*delz;
end

end