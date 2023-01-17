function vplPlot(prediction)
% vplPlot
% 
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
    ylab = 'Prediction profile likelihood';
else
    x = ar.vpl.results.z;
    y = ar.vpl.results.chi2;
    ylab = 'Validation profile likelihood';
end

alpha = ar.vpl.config.alpha;
thresh = icdf('chi2',1-alpha,1);

h1 = plot(x,y,'LineWidth',1);
hold on
line([min(x),max(x)],[thresh,thresh],'Color','k',...
    'LineWidth',1,'LineStyle','--');
xlim([min(x),max(x)]);
xlabel(ar.vpl.name,...
    'Interpreter','None');
text(mean([min(x),max(x)]),thresh,[num2str(100*(1-alpha)),'% threshold'],...
    'VerticalAlignment','Bottom','HorizontalAlignment','Center')

yl = ylim;
if nanmin(y)>-1e-5 % try to use ylim(1)=0
    yl(1)=0;
    ylim(yl)
end

if prediction==0
    h = patch([ar.vpl.results.VCI(1),ar.vpl.results.VCI(2),ar.vpl.results.VCI(2),ar.vpl.results.VCI(1)],...
        [yl(1),yl(1),yl(2),yl(2)],zeros(1,4),'FaceColor',zeros(1,3),'EdgeColor','none','FaceAlpha',0.2);

elseif prediction==1
    h = patch([ar.vpl.results.PCI(1),ar.vpl.results.PCI(2),ar.vpl.results.PCI(2),ar.vpl.results.PCI(1)],...
        [yl(1),yl(1),yl(2),yl(2)],zeros(1,4),'FaceColor',zeros(1,3),'EdgeColor','none','FaceAlpha',0.2);

elseif prediction == 2 %if both profiles should be plotted
   h2 = plot(ar.vpl.results.pred,ar.vpl.results.ppl,'LineWidth',1);
   ylab = 'Profile Value';
   legend([h1,h2],{'VPL','PPL'});
end

ylabel(ylab);
set(gca,'FontSize',12)
grid on
hold off

end

