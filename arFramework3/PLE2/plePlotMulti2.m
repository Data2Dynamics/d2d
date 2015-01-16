% Plot compact profile Likelihoods
%
% plePlotMult2i(jks, savetofile)
%
% indices           plot only parameters jks                        [all]
% savetofile        save plot                                       [false]

function plePlotMulti2(jks, savetofile)

global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end
if(~exist('jks','var') || isempty(jks))
    jks = 1:length(pleGlobals.ps);
end
if(~exist('savetofile','var'))
    savetofile = false;
end

sumples = 0;
for j=jks
    if(~isempty(pleGlobals.ps{j}))
        sumples = sumples + 1;
    end
end

if(pleGlobals.plot_point && ~pleGlobals.plot_simu)
    strtitle = sprintf('profile log-likelihood (point-wise)');
elseif(~pleGlobals.plot_point && pleGlobals.plot_simu)
    strtitle = sprintf('profile log-likelihood (simultaneous)');
else
    strtitle = sprintf('profile log-likelihood');
end

h = arRaiseFigure(pleGlobals, 'fighandel_multi2', strtitle);
set(h, 'Color', [1 1 1]);

count = 1;

minchi2 = Inf;
for jk=jks
    if(~isempty(pleGlobals.ps{jk}))
        minchi2 = min([minchi2 min(pleGlobals.chi2s{jk})]);
    end
end

% limits
dchi2 = pleGlobals.dchi2_point;
if(pleGlobals.plot_simu)
    dchi2 = pleGlobals.dchi2;
end

labels = {};

cs = lines(sumples);
for jk=jks
    if(~isempty(pleGlobals.ps{jk}))
        
        ps = pleGlobals.ps{jk};
        chi2s = pleGlobals.chi2s{jk};
        
        % profile
        plot(ps(:,jk), chi2s, 'Color', cs(count,:))
        hold on
        
        count = count + 1;
        
        labels{end+1} = arNameTrafo(pleGlobals.p_labels{jk}); %#ok<AGROW>
    end
end

% thresholds
if(pleGlobals.plot_point)
    plot(xlim, [0 0]+minchi2+chi2inv(1-pleGlobals.alpha_level, 1), 'k--')
end
if(pleGlobals.plot_simu)
    plot(xlim, [0 0]+minchi2+chi2inv(1-pleGlobals.alpha_level, pleGlobals.dof), 'k--')
end

if(count == 1)
    if(pleGlobals.plot_point && ~pleGlobals.plot_simu)
        labels{end+1} = sprintf('%2i%% (point-wise)', (1-pleGlobals.alpha_level)*100);
    elseif(~pleGlobals.plot_point && pleGlobals.plot_simu)
        labels{end+1} = sprintf('%2i%% (simultaneous)', (1-pleGlobals.alpha_level)*100);
    else
        labels{end+1} = sprintf('%2i%% (point-wise)', (1-pleGlobals.alpha_level)*100);
        labels{end+1} = sprintf('%2i%% (simultaneous)', (1-pleGlobals.alpha_level)*100);
    end
end

arSubplotStyle(gca);

ylimmax = pleGlobals.chi2+dchi2*1.3;
ylim([minchi2-dchi2*0.1 ylimmax]);

xlabel('log_{10}(parameter value)')
ylabel(pleGlobals.ylabel)

legend(labels, 'Location','EastOutside');

hold off

% save
if(savetofile && exist(pleGlobals.savePath, 'dir'))
    pleGlobals.figPathMulti2{jk} = [pleGlobals.savePath '/multi_plot2'];
    saveas(gcf, [pleGlobals.savePath '/multi_plot2'], 'fig')
    print('-depsc2', [pleGlobals.savePath '/multi_plot2']);
end



