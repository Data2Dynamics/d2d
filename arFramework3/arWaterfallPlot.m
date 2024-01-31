% arWaterfallPlot
% 
%   Basic waterfall plot, i.e. plotting sorted chi2 of the individual runs
%   of multistart optimization 
% 
%   A stair-like plot indicates that optimization works
% 
%   It can also be checked, how often the global minimum has been found.
% 
%   The term "waterfall plot" originates from the modelling workshop
%   "ODE Modelling in Systems Biology - Current Methodology and Solutions
%   of Application Issues - " that took place in Freiburg in Sep, 2017
% 
% see also arPlotChi2s

function [fig, axis] = arWaterfallPlot
global ar

fig = figure();
if(min(ar.chi2s)>0)
    axis = semilogy(1:length(ar.chi2s),sort(ar.chi2s),'.','MarkerSize',12);
else
    axis = plot(1:length(ar.chi2s),sort(ar.chi2s),'.','MarkerSize',12);
end    
xlabel('Fit index [sorted]')
ylabel('Merit after optimization')
xlim([0,length(ar.chi2s)+1])
title(['Waterfall plot: ',num2str(length(ar.chi2s)),' fits, ',num2str(sum(isinf(ar.chi2s))), ' failed'])
grid on
