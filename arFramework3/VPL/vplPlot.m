function vplPlot(prediction)
% vplPlot(prediction)
%
% Basic visualization of the calculated validation profile.
%
% prediction: Plot Prediction instead of validation profile [false]
%
% The default for showing the prediction/validation profile is controlled
% by the same flag which decides whether prediction or validation profile
% has been calculated (ar.vpl.config.prediction).
%
% See also: InitVPL, VPL 

global ar

if ~exist('prediction','var')
    prediction = ar.vpl.config.prediction;
end

if (prediction)
    x = ar.vpl.results.pred;
    y = ar.vpl.results.ppl;
    ylab = 'PPL';
else
    x = ar.vpl.results.z;
    y = ar.vpl.results.chi2;
    ylab = 'VPL';
end

thresh = icdf('chi2',0.95,1);

figure
hold on
plot(x,y,'LineWidth',1);
line([min(x),max(x)],[thresh,thresh],'Color','k',...
    'LineWidth',1,'LineStyle','--');
xlim([min(x),max(x)]);
xlabel(ar.model(ar.vpl.general.m).data(ar.vpl.general.d).y(ar.vpl.general.idpred),...
    'Interpreter','None');
ylabel(ylab);
text(mean([min(x),max(x)]),thresh,'95% treshold','VerticalAlignment','Bottom',...
    'HorizontalAlignment','Center')
hold off

end

