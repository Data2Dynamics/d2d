function vplPlot(prediction)
% vplPlot(prediction)
%
% Basic visualization of the calculated validation profile.
%
% prediction: 0 - Plot Validation Profile
%             1 - Plot Prediction Profile
%             2 - Plot both
%               [default: ar.vpl.config.prediction]
%
% The default for showing the prediction/validation profile is controlled
% by the same flag which decides whether prediction or validation profile
% has been calculated (ar.vpl.config.prediction).
%
% See also: InitVPL, VPL 

global ar

if ~exist('prediction','var')
    % Plot what was intended to be calculated:
    prediction = ar.vpl.config.prediction;
end

if prediction == 1
    x = ar.vpl.results.pred;
    y = ar.vpl.results.ppl;
    ylab = 'PPL';
else
    x = ar.vpl.results.z;
    y = ar.vpl.results.chi2;
    ylab = 'VPL';
end

alpha = ar.vpl.config.alpha;
thresh = icdf('chi2',1-alpha,1);

figure
hold on

h1 = plot(x,y,'LineWidth',1);
line([min(x),max(x)],[thresh,thresh],'Color','k',...
    'LineWidth',1,'LineStyle','--');
xlim([min(x),max(x)]);
xlabel(ar.model(ar.vpl.general.m).data(ar.vpl.general.d).y(ar.vpl.general.idpred),...
    'Interpreter','None');
text(mean([min(x),max(x)]),thresh,[num2str(100*(1-alpha)),'% threshold'],...
    'VerticalAlignment','Bottom','HorizontalAlignment','Center')

if prediction == 2 %if both profiles should be plotted
   h2 = plot(ar.vpl.results.pred,ar.vpl.results.ppl,'LineWidth',1);
   ylab = 'Profile Value';
   legend([h1,h2],{'VPL','PPL'});
end

ylabel(ylab);

hold off

end

