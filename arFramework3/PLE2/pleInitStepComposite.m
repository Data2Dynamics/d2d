% adaptive step reconciling chi^2 increase
% direct & linear step methods

function [pStep, dpNew, beta, alpha] = pleInitStepComposite(jk, pLast, dpLast)

[pStep1, dpNew1, beta1, alpha1] = pleInitStepDirect(jk, pLast, dpLast);
[pStep2, dpNew2, beta2, alpha2] = pleInitStepLinear(jk, pLast, dpLast);

index = 1;
if(index==1)
    pStep = pStep1;
    dpNew = dpNew1;
    beta = beta1;
    alpha = alpha1;
else
    pStep = pStep2;
    dpNew = dpNew2;
    beta = beta2;
    alpha = alpha2;
end
