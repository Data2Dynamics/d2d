% arPlotParameterHists([ps], [jks], [nbins])
%
% Plot histograms of parameter matrices (a row corresponds to one parameter set)
% 
%   ps     [ar.ps]                               Parameter matrix
%   jks    [find(ar.qDynamic==1 & ar.qFit==1)]   Only ps(:,jks) is plotted
%   nbins  [50]                                  Number of bins of the histogram
%
%
% Example1:
% arFitLHS(100)
% arPlotParameterHists(ar.ps)
% 
% Example2: Plotting of parameters with failed integrations.
% arChi2LHS(100)
% arPlotParameterHists(ar.ps_errors)
%
% See also: arPlotParameterPatterns and arPlotParameterProfiles

function arPlotParameterHists(ps, jks, nbins)

global ar

if(~exist('ps','var'))
    ps = ar.ps;
end

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

if(~exist('nbins','var'))
    nbins = 50;
end

figure(1); clf;

[nrows, ncols] = arNtoColsAndRows(length(jks));

count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    arSubplotStyle(g)
    hist(ps(:,j), nbins);
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r');
    hold off
    if mod(count,ncols) == 1
        ylabel('count')
    end
    if count > (nrows-1) * ncols
        xlabel('log10(p)')
    end
    count = count + 1;
    title(arNameTrafo(ar.pLabel{j}));
end







