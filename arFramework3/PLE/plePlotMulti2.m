% plePlotMult2i([jks], [savetofile])
%
% Plot multiple profile Likelihoods in subplots
%
% jks               parameter indices to be plotted                 [1:length(ar.ple.ps)]
% savetofile        save plot                                       [false]
% 
% See also arPLEInit, ple, plePlot, arPlotMulti

function plePlotMulti2(jks, savetofile)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end
if(~exist('jks','var') || isempty(jks))
    jks = 1:length(ar.ple.ps);
end
if(~exist('savetofile','var'))
    savetofile = false;
end

sumples = 0;
for j=jks
    if(~isempty(ar.ple.ps{j}))
        sumples = sumples + 1;
    end
end

if(ar.ple.plot_point && ~ar.ple.plot_simu)
    strtitle = sprintf('profile log-likelihood (point-wise)');
elseif(~ar.ple.plot_point && ar.ple.plot_simu)
    strtitle = sprintf('profile log-likelihood (simultaneous)');
else
    strtitle = sprintf('profile log-likelihood');
end

h = myRaiseFigure(strtitle);
set(h, 'Color', [1 1 1]);

count = 1;

minchi2 = Inf;
for jk=jks
    if(~isempty(ar.ple.ps{jk}))
        minchi2 = min([minchi2 min(ar.ple.chi2s{jk})]);
    end
end

% limits
dchi2 = ar.ple.dchi2_point;
if(ar.ple.plot_simu)
    dchi2 = ar.ple.dchi2;
end

labels = {};

cs = lines(sumples);
for jk=jks
    if(~isempty(ar.ple.ps{jk}))
        
        ps = ar.ple.ps{jk};
        chi2s = ar.ple.chi2s{jk};
        
        % profile
        plot(ps(:,jk), chi2s, 'Color', cs(count,:))
        hold on
        
        count = count + 1;
        
        labels{end+1} = arNameTrafo(ar.ple.p_labels{jk}); %#ok<AGROW>
    end
end

% thresholds
if(ar.ple.plot_point)
    plot(xlim, [0 0]+minchi2+arChi2inv(1-ar.ple.alpha_level, 1), 'k--')
end
if(ar.ple.plot_simu)
    plot(xlim, [0 0]+minchi2+arChi2inv(1-ar.ple.alpha_level, ar.ple.dof), 'k--')
end

if(count == 1)
    if(ar.ple.plot_point && ~ar.ple.plot_simu)
        labels{end+1} = sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100);
    elseif(~ar.ple.plot_point && ar.ple.plot_simu)
        labels{end+1} = sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100);
    else
        labels{end+1} = sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100);
        labels{end+1} = sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100);
    end
end

arSubplotStyle(gca);

ylimmax = ar.ple.merit+dchi2*1.3;
ylim([minchi2-dchi2*0.1 ylimmax]);

xlabel('log_{10}(parameter value)')
ylabel(ar.ple.ylabel)

legend(labels, 'Location','EastOutside');

hold off

% save
if(savetofile && exist(ar.ple.savePath, 'dir'))
    ar.ple.figPathMulti2{jk} = [ar.ple.savePath '/multi_plot2'];
    saveas(gcf, [ar.ple.savePath '/multi_plot2'], 'fig')
    print('-depsc2', [ar.ple.savePath '/multi_plot2']);
end


function h = myRaiseFigure(figname)

openfigs = get(0,'Children');

figcolor = [1 1 1];

if(isfield(ar.ple, 'fighandel_multi2') && ~isempty(ar.ple.fighandel_multi2) && ...
    ar.ple.fighandel_multi2 ~= 0 && ...
    sum(ar.ple.fighandel_multi2==openfigs)>0 && ...
    strcmp(get(ar.ple.fighandel_multi2, 'Name'), figname))

    h = ar.ple.fighandel_multi2;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.1 0.6 0.8]);
    set(h,'Color', figcolor);
    ar.ple.fighandel_multi2 = h;
end


