% This function plots the fit results
function plotFitTransient(fit,folder)
if ~exist('folder','var') || isempty(folder)
    folder = 'plotFitTransient';
end

if iscell(fit)
    if ~isdir(folder)
        mkdir(folder);
    end
    
    for i=1:length(fit)
        try
            plotFitTransient(fit{i});
            print(sprintf('%s%s%s',folder,filesep,fit{i}.label),'-dpng');
        catch ERR
            warning('Error for %s ',fit{i}.label);
            warning(ERR.message)
        end
    end
else       
    clf('reset')
    if ~isfield(fit.data,'ystd')
    patch([fit.tFine;fit.tFine(end:-1:1)],...
        [fit.yFineSimu+fit.ystdFineSimu;fit.yFineSimu(end:-1:1)-fit.ystdFineSimu(end:-1:1)],...
        zeros(length(fit.tFine)*2,1)-2*eps,'FaceColor',0.7*ones(1,3),'EdgeColor',0.7*ones(1,3));
    hold on
    plot(fit.data.tExp,fit.data.yExp,'ko','LineWidth',2)
    else
        hold on
        errorbar(fit.data.tExp,fit.data.yExp,fit.data.ystd,'ko')
    end
    plot(fit.tFine,fit.yFineSimu,'-','LineWidth',1,'Color',0*ones(1,3));
    set(gca,'FontSize',14,'LineWidth',2)
    axis tight
        
    xl = xlim;
    xl = xl + diff(xl)*[-.03,.03];
    xlim(xl);
    if range(fit.data.yExp)>2e-6
        yl = ylim;
        yl = yl + diff(yl)*[-.03,.03];
        ylim(yl);
    else
        ylim([-.1,1.1])
    end
    title(strrep(fit.label,'_','\_'));
    saveas(gcf,sprintf('%s',folder),'png'); % saveas overwrites image
end
