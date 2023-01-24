% plePlotNice(whichOne,figname)
% 
%   Nice graphs of individual profiles
% 
% Examples:
% plePlotNice('KD')
% plePlotNice('KD','ModelName') % with title and saving as png
% plePlotNice({'KD','hill'}) % several params
% plePlotNice(1:3)      % several params
% plePlotNice           % all params


function plePlotNice(whichOne,figname)
global ar
if ~exist('whichOne','var') || isempty(whichOne)
    whichOne = 1:length(ar.ple.ps);
end
if ~exist('figname','var') || isempty(figname)
    figname = '';
end

if ischar(whichOne) || iscell(whichOne)
    [~,ind] = intersect(ar.pLabel,whichOne,'stable');
else
    ind = whichOne;
end

if length(ind)>1 % many single profiles:
    for i=1:length(ind)
        plePlotNice(ind(i),figname)
    end
else % single profile:
    if isempty(ar.ple.chi2s{ind})
        return  % profile not calculated
    end
    chi2s = ar.ple.chi2s{ind};
    ps = ar.ple.ps{ind}(:,ind);
    notnan = find(~isnan(chi2s));
    chi2s = chi2s(notnan);
    ps = ps(notnan);
    chi2min = nanmin(chi2s);
    
    thresh = arChi2inv(1-ar.ple.alpha_level, 1);
    
    imin = find(chi2s==chi2min);
    imin = round(mean(imin)); % if not unique use the middle
    chi2left = chi2s(1:imin);
    plb = interp1(chi2left,ps(1:imin),chi2min+thresh);
    chi2right = chi2s(imin:end);
    pub = interp1(chi2right,ps(imin:end),chi2min+thresh);
    
    figure
    plot(ps,chi2s,'k','LineWidth',2)
    set(gca,'FontSize',14)
    grid on
    if ar.ple.qLog10(ind)==1
        xlabel([strrep(ar.pLabel{ind},'_','\_'),'  [log10]'])
    else
        xlabel(strrep(ar.pLabel{ind},'_','\_'))
    end
    ylabel('-2 log likelihood')
    text(mean(xlim), nanmin(chi2s)+thresh,...
        sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100), 'Color', 'k', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
    hold on
    xl = xlim;
    xlim(xl);    
    yl = ylim;
    yl(1) = chi2min;
    ylim(yl);
    h = patch([plb,pub,pub,plb],[yl(1),yl(1),yl(2),yl(2)],zeros(1,4),'FaceColor',zeros(1,3),'EdgeColor','none','FaceAlpha',0.3);
    plot(xl,nanmin(chi2s)+arChi2inv(1-ar.ple.alpha_level, 1)*ones(1,2),'k--','LineWidth',2);
    if ar.ple.qLog10(ind)
        title(sprintf('Estimate=%g, CI=[%g, %g]',10^ar.ple.p(ind),10^ar.ple.conf_lb_point(ind),10^ar.ple.conf_ub_point(ind)))
    else
        title(sprintf('Estimate=%g, CI=[%g, %g]',ar.ple.p(ind),ar.ple.conf_lb_point(ind),ar.ple.conf_ub_point(ind)))
    end
    if ~isempty(figname)
%         title(str2label(figname))
        print([figname,'_',ar.pLabel{ind},'_PL'],'-dpng')
    end
end

