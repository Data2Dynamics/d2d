% arLogPlotting(logplotting) 
% 
% Set the log-plotting option of all plots 
% 
%  logplotting     Boolean, false: plot on linear scale, 
%                           true:  plot on log-scale
%

function arLogPlotting(logplotting)
global ar

for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        for jd = 1:length(ar.model(jm).data)
            ar.model(jm).data(jd).logplotting(:) = logplotting;
        end
    end
end