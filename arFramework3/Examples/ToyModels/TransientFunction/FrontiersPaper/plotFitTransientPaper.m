% plotFitTransientPaper(fit,folder,fitODE)
% 
% This function plots the fit results in the paper format
% 
%   fit     struct or cell of structs as generated containing the RTF
%           result
% 
%   fitODE  result for ODE model
% 
%   See Setup_x in the paper analysis (Bachmann folder) or
%   Analysis_CollectFits.m for examples
% 
% Examples:
% load result_Paper
% plotFitTransientPaper(res,'plotFitTransientPaper_Data',sim)
% plotFitTransientPaper(res,'plotFitTransientPaper_Data_woODE')

function plotFitTransientPaper(fit,folder,fitODE)
if ~exist('folder','var') || isempty(folder)
    folder = 'plotFitTransientPaper';
end
if ~exist('fitODE','var') || isempty(fitODE)
    if iscell(fit)
        fitODE = cell(size(fit));
    else
        fitODE = [];
    end
end


if iscell(fit)
    if ~isfolder(folder)
        mkdir(folder);
    end
    
%     for i=[14,22,33:35,42,44,47,54,75,78,101,103:104,254,277,300,305,323,339:346,362:368,385:391,401,408:414]%1:length(fits)    
    for i=1:length(fit)
%         try
            plotFitTransientPaper(fit{i},folder,fitODE{i});  
    
%             set(gcf,'PaperSize',[10,10])
%             set(gcf,'PaperPosition', [0 0 10 10]);
%             print(sprintf('%s%s%3i_%s',folder,filesep,i,fit{i}.label),'-dpdf');
            print(sprintf('%s%s%3i_%s',folder,filesep,i,fit{i}.label),'-dpng');

%         catch ERR
%             warning('Error for %s ',fit{i}.label);
%             warning(ERR.message)
%         end
    end
else       
    
    h4 = [];
    
    clf('reset')
    L = lines;
    if ~isfield(fit.data,'ystd')
        h3=patch([fit.tFine;fit.tFine(end:-1:1)],...
            [fit.yFineSimu+fit.ystdFineSimu;fit.yFineSimu(end:-1:1)-fit.ystdFineSimu(end:-1:1)],...
            zeros(length(fit.tFine)*2,1)-2*eps,'FaceColor',0.7*ones(1,3),'EdgeColor',0.7*ones(1,3));
        hold on
        h1=plot(fit.data.tExp,fit.data.yExp,'k-','LineWidth',2);
    else
        hold on
        h1 = errorbar2(fit.data.tExp,fit.data.yExp,fit.data.ystd,'ko');
        h3 = plot(fit.data.tFine,fit.data.yFine,'k','LineWidth',2);
    end
    if ~isempty(fitODE)
        h4=plot(fitODE.tFine,fitODE.yFineFit,'-','LineWidth',2,'Color',L(2,:));
    end
    h2=plot(fit.tFine,fit.yFineSimu,'-','LineWidth',2,'Color',L(1,:));
        
    set(gca,'FontSize',14,'LineWidth',2)
    axis tight
    xlabel('time  [min]')
    ylabel('concentration  [a.u.]')
        
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
    tit = fit.label;
    title(tit);

    yl = ylim;
    
    y = fit.data.yExp;
    y = y(ceil(length(y)*0.7):end);
    if max(y)>yl(1)+0.7*range(yl)   
        if isempty(h4) && isfield(fit,'rmse')
            legend([h2(1)],sprintf('RTF (RMSE=%.2f)',fit.rmse),'Location','SouthEast');
        elseif ~isfield(fit.data,'ystd')
           if isfield(fit,'rmse')
               legend([h1(1),h2(1),h3(1)],sprintf('ODE (RMSE=%.2f)',fitODE.rmse),sprintf('RTF (RMSE=%.2f)',fit.rmse),'Approx. error','Location','SouthEast');
           else
               legend([h1(1),h2(1),h3(1)],'ODE','RTF',sprintf('Approx. error = %.2f',fit.approxErr),'Location','SouthEast');
           end
        else
            if isfield(fit,'rmse')
                legend([h3(1),h1(2),h2(1),h4(1)],'True dynamics','Data',sprintf('RTF fit (RMSE=%.2f)',fit.rmse),sprintf('ODE fit (RMSE=%.2f)',fitODE.rmse),'Location','SouthEast');
            else
                legend([h3(1),h1(2),h2(1),h4(1)],'True dynamics','Data','RTF','ODE fit','Location','SouthEast');
            end
        end
    else
        if isempty(h4) && isfield(fit,'rmse')
            legend([h2(1)],sprintf('RTF (RMSE=%.2f)',fit.rmse));
        elseif ~isfield(fit.data,'ystd') 
            legend([h1(1),h2(1),h3(1)],'ODE','RTF',sprintf('Approx. error = %.2f',fit.approxErr));
        else
            if isfield(fit,'rmse')
                legend([h3(1),h1(2),h2(1),h4(1)],'True dynamics','Data',sprintf('RTF fit (RMSE=%.2f)',fit.rmse),sprintf('ODE fit (RMSE=%.2f)',fitODE.rmse));
            else
                legend([h3(1),h1(2),h2(1),h4(1)],'True dynamics','Data','RTF','ODE fit');
            end
            
        end
    end
end
