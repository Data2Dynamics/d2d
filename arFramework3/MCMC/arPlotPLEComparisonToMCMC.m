% arPlotPLEComparisonToMCMC([ParaPerFigure,ParaVector,BayesianCredibleInterval,LogScale])
%
%  ParaPerFigure [16]
%       Number of parameters that are plotted into same figure 
%
%  ParaVector [find(ar.qFit==1)]
%       Give parameter vector to only plot correlations of certain
%       parameters, e.g. ParaVector = [ 1 2 3 4 6 9 112 114]
%
%  BayesianCredibleInterval (true) 
%       Specify, if bayesian 95% credible interval is plotted
%
%  LogScale (true)
%       Specify, whether marginalized (histogram) is plotted in log scale
%       to accomodate for the fact that it is proportional to the
%       likelihood, not the log-Likelihood



function arPlotPLEComparisonToMCMC(varargin)

global ar


if(isempty(ar.ple.ps))
    error('Please run arPLEInit and ple first!, No PLE data found!')
end
    

if nargin > 3
    LogScale = varargin{4};    
else 
    LogScale = true;
end
if nargin > 2
    BayesianCredibleInterval = varargin{3};    
else 
    BayesianCredibleInterval = true;
end


if nargin > 1 && (~isempty(varargin{2}))
    ParaVector = varargin{2};
else 
    ParaVector = find(ar.qFit==1);
end

if nargin > 0
    ParaPerFigure = varargin{1};    
else
    ParaPerFigure = 16;
end



ConfLevel = 0.05;

conf_lb_point = ar.ple.conf_lb_point;
conf_ub_point = ar.ple.conf_ub_point;
pBestFit = ar.ple.p;

RunNumber = size(ar.ps,1);

NumberOfFigures = ceil(length(ParaVector)/ParaPerFigure);

for FigureNumber = 1:NumberOfFigures

    if FigureNumber < NumberOfFigures    
        FigureParaVector = ParaVector(1+(FigureNumber-1)*ParaPerFigure:(FigureNumber-1)*ParaPerFigure+ParaPerFigure);
    else
        FigureParaVector = ParaVector(1+(FigureNumber-1)*ParaPerFigure:end);
    end
    if length(FigureParaVector) > 8
        qq = ceil(length(FigureParaVector)/2);
    else
        qq = length(FigureParaVector);
    end
    NumberOfColumns = min([qq 8]);
    NumberOfRows = ceil(2*length(FigureParaVector)/NumberOfColumns);

    figure('Name',sprintf('PLE compared to Marginalized - Figure  %d / %d',FigureNumber,NumberOfFigures),'NumberTitle','off')
    for paraindex=FigureParaVector
        hold on
        i = find(paraindex==FigureParaVector);
          
        subplot(NumberOfRows,NumberOfColumns,i+floor((i-1)/NumberOfColumns)*NumberOfColumns)       % add first plot in 2 x 1 grid

        %h1 = histogram(ar.ps(:,paraindex));
        h1 = histogram(ar.ps(:,paraindex),100);
        set(h1, 'FaceColor','none');
        
        %         % Log Scale of Histogram to make it comparable to PLE
        if LogScale == true
            set(gca,'YScale','log')
            limsy=get(gca,'YLim');
            set(gca,'Ylim',[.7 limsy(2)]);     
        end
        
        if ar.qLog10(paraindex) == 1
            xlabel(['log_{10}(' arNameTrafo(ar.pLabel{paraindex}) ')'])
        else
            xlabel([ arNameTrafo(ar.pLabel{paraindex})])
        end
        ylabel('posterior PDF');
        
        % Bayesian confidence interval 95 %
        if BayesianCredibleInterval == true
            Sorted = sort(ar.ps(:,paraindex));
            LowerConfidenceIndex = round(RunNumber*ConfLevel/2);
            UpperConfidenceIndex = round(RunNumber-RunNumber*ConfLevel/2);
            line([Sorted(LowerConfidenceIndex), Sorted(LowerConfidenceIndex)], ylim, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
            line([Sorted(UpperConfidenceIndex), Sorted(UpperConfidenceIndex)], ylim, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');          
        end

        title(ar.pLabel{paraindex});
        ax = gca; % current axes


        hold on
        subplot(NumberOfRows,NumberOfColumns,NumberOfColumns+i+floor((i-1)/NumberOfColumns)*NumberOfColumns)


        show_hit_bound = 1:length(ar.ple.ps);
        plot_hit_bound = true;
        plot_thresholds = true;
        minchi2 = Inf;
        if(~isempty(ar.ple.ps{paraindex}))
            minchi2 = min([minchi2 min(ar.ple.chi2s{paraindex})]);
        end

        ps = ar.ple.ps{paraindex};
        chi2s = ar.ple.chi2s{paraindex};

        qCloseToUB = ps > ones(length(chi2s),1) * (ar.ub - ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * ar.qFit==1;
        qCloseToLB = ps < ones(length(chi2s),1) * (ar.lb + ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * ar.qFit==1;

        qhitbound = false(size(ps));
        qhitbound(:,ar.qFit==1) = ar.ple.gradient{paraindex}(:,ar.qFit==1) > ar.ple.grad_thres & qCloseToLB(:,ar.qFit==1) | ...
            ar.ple.gradient{paraindex}(:,ar.qFit==1) < -ar.ple.grad_thres & qCloseToUB(:,ar.qFit==1);

        % profile
        qplot = true(size(ps,1),1);
        if(sum(show_hit_bound==paraindex)>0)
            plot(ps(:,paraindex), chi2s, 'k-', 'LineWidth', 1)
        else
            qplot = sum(qhitbound,2)==0;
            plot(ps(qplot,paraindex), chi2s(qplot), 'k-', 'LineWidth', 1)
        end
        hold on
        % boundary values
        if(sum(show_hit_bound==paraindex)>0 && plot_hit_bound)
            plot(ps(sum(qhitbound,2)>0,paraindex), chi2s(sum(qhitbound,2)>0), 'ko', 'LineWidth', 1)
        end

        % limits
        dchi2 = ar.ple.dchi2_point;
        if(ar.ple.plot_simu)
            dchi2 = ar.ple.dchi2;
        end

        ylimmax = ar.ple.merit+dchi2*1.3;
        if(plot_thresholds)
            ylim([minchi2-dchi2*0.1 ylimmax]);
        end

        qbelowchi2 = chi2s < ylimmax;
        xlimtmp2 = (max(ps(qbelowchi2,paraindex))-min(ps(qbelowchi2,paraindex)))*0.05;
        if(xlimtmp2>0)
            if(show_hit_bound)
                xlimtmp = [min(ps(qbelowchi2,paraindex))-xlimtmp2 max(ps(qbelowchi2,paraindex))+xlimtmp2];
            else
                xlimtmp = [min(ps(sum(qhitbound,2)==0 & qbelowchi2,paraindex))- ...
                    xlimtmp2 max(ps(sum(qhitbound,2)==0 & qbelowchi2,paraindex))+xlimtmp2];
            end
            xlim(xlimtmp);
        end
        
        
        if ar.qLog10(paraindex) == 1
            xlabel(['log_{10}(' arNameTrafo(ar.pLabel{paraindex}) ')'])
        else
            xlabel([ arNameTrafo(ar.pLabel{paraindex})])
        end
        ylabel(ar.ple.ylabel)



        % thresholds

        if(plot_thresholds)

            if(ar.ple.plot_point)
                plot(xlim, [0 0]+minchi2+chi2inv(1-ar.ple.alpha_level, 1), 'r-')
            end
            if(ar.ple.plot_simu)
                plot(xlim, [0 0]+minchi2+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), 'r-')
            end



            line([conf_lb_point(paraindex), conf_lb_point(paraindex)], ylim, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'r');
            line([conf_ub_point(paraindex), conf_ub_point(paraindex)], ylim, 'LineStyle', '--', 'LineWidth', 1, 'Color', 'r');

        end
        %%%%%%%%%%%%%% END OF PLOTTING RECALCULATED
        %%%%%%%%%%%%%% THRESHOLDS

        bx = gca;
        bx.YDir = 'reverse';
        xlim(ax.XLim);            
    end   
    set(gcf, 'Position', get(0, 'Screensize'));
end

end