% RMSE after scaling the range to 1.

function [rmse,ytrue2,yest2] = CalculateRMSE(ytrue,yest,SD)

m = min(ytrue);
if exist('SD','var') && ~isempty(SD)
    fac = SD;
else
    fac = range(ytrue);
end

ytrue2 = (ytrue-m)/fac;
yest2 = (yest-m)/fac;

rmse = sqrt(mean((ytrue2 - yest2).^2));
