% RMSE of structs as calculated by Setup_x
%   ytrue,yest  scaled to range equals to one

function [rmse,ytrue,yest] = CalculateRmseSim(sim)

rmse = NaN(size(sim));
ytrue = NaN(length(sim{1}.yFine),length(rmse));
yest = NaN(length(sim{1}.yFine),length(rmse));
for i=1:length(sim)
    SD = mean(sim{i}.ystd);
% try        
    [rmse(i),ytrue(:,i),yest(:,i)] = CalculateRMSE(sim{i}.yFine,sim{i}.yFineFit,SD);    
% catch
%     i
% end
end
