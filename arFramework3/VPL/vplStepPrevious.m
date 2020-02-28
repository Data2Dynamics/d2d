function delz = vplStepPrevious(z_old,chi2dif,delz,gen_struct)
% delz = vplStepPrevious(z_old,chi2dif,delz,gen_struct)
%
% Does the validation data step. Primitive version, no dynamic step size
% control, but adaption based on the last step.

%Magic factors:
maxstepsize       = gen_struct.temp.maxstepsize;
minstepsize       = gen_struct.temp.minstepsize;
stepfactor        = gen_struct.temp.stepfactor;
chi2dif_max       = gen_struct.temp.chi2dif_max;
chi2dif_min       = gen_struct.temp.chi2dif_min;

%Adjust delz based on last step step:
if abs(chi2dif) > chi2dif_max
    if delz/stepfactor > minstepsize
        delz = delz/stepfactor;
    end
elseif abs(chi2dif) < chi2dif_min
    if delz*stepfactor < maxstepsize
        delz = delz*stepfactor;
    end
end

end