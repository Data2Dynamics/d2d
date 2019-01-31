% arSetLogPlotting(logplotting, jm, jp) 
% 
% Set the log-plotting option of the plot in ar.model(jm).plot(jp)
% 
%   logplotting        Boolean, false: plot on linear scale, 
%                               true:  plot on log-scale
%   jm          [1]    Index of the model
%   jp          [1]    Index of the plot


function arSetLogPlotting(logplotting, jm, jp) 
global ar

if(~exist('jm','var') || isempty(jm) )
    jm = 1;
end

if(~exist('jp','var') || isempty(jp) )
    jp = 1;
end

for jd = ar.model(jm).plot(jp).dLink
    ar.model(jm).data(jd).logplotting(:) = logplotting;
end
