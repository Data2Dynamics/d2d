% arPlotParameterProfiles([jks])
%
% Plots the values of ar.ps(:,jks) versus ar.chi2s(:,jks) relative to the
% best fit.
%
%   jks         [find(ar.qDynamic==1 & ar.qFit==1)]      Index of parameters to be plotted
%
%
%  To analyse the results of arFitLHS or arChi2LHS
%
% Example:
% arFitLHS(100)
% arPlotParameterProfiles([1 2 4])  % plot distribution of parameters ar.p([1 2 4])
%
% See also: arPlotParameterPatterns and arPlotParameterHists


function arPlotParameterProfiles(jks)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qDynamic==1 & ar.qFit==1);
end

figure(1); clf;

dchi2 = chi2inv(0.95, 1);

[nrows, ncols] = arNtoColsAndRows(length(jks));
ncols
count = 1;
for j=jks
    g = subplot(nrows,ncols,count);
    arSubplotStyle(g)

    plot(ar.ps(:,j), log10(ar.chi2s-min(ar.chi2s)+1), 'xk');
    xlim([ar.lb(j) ar.ub(j)]);
    hold on
    plot([ar.p(j) ar.p(j)], ylim, 'r--');
    plot(xlim, [0 0], 'k--');
    plot(xlim, log10([dchi2 dchi2]), 'k:');
    hold off
    title(arNameTrafo(ar.pLabel{j}));
    if mod(count,ncols) == 1
        ylabel('log10(\Deltachi2)')
    end
    if count > (nrows-1) * ncols
        xlabel('log10(p)')
    end
    count = count + 1;
    arSpacedAxisLimits;
end


