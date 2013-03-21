% toggle log-plotting option
%
% arLogPlotting(logplotting)

function arLogPlotting(logplotting)
global ar

for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        for jd = 1:length(ar.model(jm).data)
            ar.model(jm).data(jd).logplotting(:) = logplotting;
        end
    end
end