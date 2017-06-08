%   plotopt = arWhichYplot(m,d)
%   plotopt = arWhichYplot(m,d,it,iy)
%
%   m   model index
%   d   data index
%   it  data time index [ default: 1:size(ar.model.data.yExpStd,1) ]
%   iy  data observable index [ default: 1:size(ar.model.data.yExpStd,2) ]
%
%
% This function evaluates
%   - ar.config.ploterrors
%   - ar.config.fiterrors
%   - ar.model.data.yExpStd
% and determines the way of plotting data.
%
% The function replaces multiple complex logical operations.
%
%
%   plotopt = 0 means plotting nothing
%   plotopt = 1 means plotting only data points
%   plotopt = 2 means plotting data and errorbars
%   plotopt = 3 means plotting data and error model
%   plotopt = 4 means plotting data and prediction bands
%   plotopt = 5 means plotting data, errorbars and error model
%
% Doku:
% https://github.com/Data2Dynamics/d2d/wiki/Plotting-options-and-the-meaning-of-ar.config.ploterrors


function plotopt = arWhichYplot(m,d,it,iy)
global ar

if ~exist('it','var') || isempty(it)
    it = [];
end
if ~exist('iy','var') || isempty(iy)
    iy = [];
end

if length(d)>1  % occurs in case of Dose-Response
    plotopt = NaN(size(d));
    for i=1:length(d)
        plotopt(i) = arWhichYplot(m,d(i),it,iy);
    end
    plotopt = unique(plotopt);
    if length(plotopt)>1
        if sum(plotopt==4)>0
            plotopt = 4; % prediction bands
        elseif sum(plotopt==5)>0
            plotopt = 5;  % errorbar and error bands
        elseif sum(plotopt==2)>0 && sum(plotopt==3)>0
            plotopt = 5;
        else
            plotopt = max(plotopt);
        end
    end
else
    if isempty(it)
        it = 1:size(ar.model(m).data(d).yExpStd,1);
    end
    if isempty(iy)
        iy = 1:size(ar.model(m).data(d).yExpStd,2);
    end
    
    if ar.config.ploterrors == -2
        plotopt = 1;
    elseif ar.config.ploterrors == -1 % prediction conf. bands
        plotopt = 4;
    elseif ar.config.ploterrors == 0  % don't externally control plotting. Instead evaluate ar.config.fiterrors and availability of exp. errors.
        if isempty(ar.model(m).data(d).yExpStd)
            plotopt = 5;
        elseif sum(isnan(ar.model(m).data(d).yExpStd(it,iy)) )==0  && ar.config.fiterrors~=1  % all exp. errors available && exp. errors used
            plotopt = 2;
        elseif (ar.config.fiterrors==1 && ar.config.ploterrors==0) || (sum(~isnan(ar.model(m).data(d).yExpStd(it,iy)) & ~isnan(ar.model(m).data(d).yExp(it,iy)))==0   && ar.config.fiterrors~=-1) % no exp. errors available && error model used
            plotopt = 3;
        else
            plotopt = 5;
        end
    elseif ar.config.ploterrors == 1  % error bars, no error model
        plotopt = 2;
    elseif ar.config.ploterrors == 2  % error model
        plotopt = 3;
    else
        error('ar.config.ploterrors = %d is not implemented.',ar.config.ploterrors);
    end
    
end