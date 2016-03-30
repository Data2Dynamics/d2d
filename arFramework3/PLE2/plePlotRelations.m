% Plot relation between parameters
%
% plePlotRelations([jks, logplot, jresponse])
%
% jks:          vector of parameter indices           []
%               if jks=[] a matrix plot is generated
% logplot:      plot logarithmically                  [false]
% jresponse:    parameter taken as response           [0]
%               if jresponse=0 jks are considered

function plePlotRelations(jks, logplot, jresponse)

global pleGlobals;

if(isempty(pleGlobals))
    error('perform ple before usage');
end
if(isempty(pleGlobals.ps))
    return
end

if(~exist('jks', 'var'))
    jks = [];
end
if(~exist('logplot', 'var'))
    logplot = false;
end
if(~exist('jresponse', 'var'))
    jresponse = 0;
end

% Setup k-Matrix
if(jresponse==0)
    k = [];
    chi2 = [];
    
    for j=find(pleGlobals.q_fit)
        if (j<=size(pleGlobals.ps,2))
            if(~isempty(pleGlobals.ps{j}))
                k = [k;pleGlobals.ps{j}]; %#ok<AGROW>
                chi2 = [chi2 pleGlobals.chi2s{j}]; %#ok<AGROW>
            end
        end
    end
else
    k = pleGlobals.ps{jresponse};
    chi2 = pleGlobals.chi2s{jresponse};
end

if(logplot)
    k = log10(k);
end

nfarben = 50;

% farben = winter(nfarben);
% farben(end,:) = [1 0 0];

farben = gray(nfarben);
farben(end,:) = [1 1 1];

farbenindex = (chi2-min(chi2))./pleGlobals.dchi2;
farbenindex = 1+floor(farbenindex*nfarben);
farbenindex(farbenindex>nfarben) = nfarben;

figure(length(pleGlobals.p)+2)

if(isempty(jks))
    plotmatrix(k(:,pleGlobals.q_fit))
elseif(length(jks) == 1)
    hist(k(:,jks), ...
        sum(~isnan(k(:,jks)))/5)
    xlabel(pleGlobals.p_labels{jks(1)})
    ylabel('frequency')
elseif(length(jks) == 2)
    subplot(1,1,1)
    for j=1:length(k)
        if(~isnan(farbenindex(j)))
            plot(k(j,jks(1)), k(j,jks(2)), 'o', ...
                'MarkerEdgeColor','k','MarkerFaceColor', farben(farbenindex(j),:), ...
                'LineWidth',0.5,'MarkerSize',4);
            hold on
        end
    end
    hold off
    xlabel(pleGlobals.p_labels{jks(1)})
    ylabel(pleGlobals.p_labels{jks(2)})
    grid on
elseif(length(jks) == 3)
    subplot(1,1,1)
    for j=1:length(k)
        if(~isnan(farbenindex(j)))
            plot3(k(j,jks(1)), k(j,jks(2)), k(j,jks(3)), '.', ...
                'MarkerEdgeColor','k','MarkerFaceColor', farben(farbenindex(j),:), ...
                'LineWidth',0.5,'MarkerSize',4);
            hold on
        end
    end
    hold off
    xlabel(pleGlobals.p_labels{jks(1)})
    ylabel(pleGlobals.p_labels{jks(2)})
    zlabel(pleGlobals.p_labels{jks(3)})
    grid on
elseif(length(jks) > 3)
    plotmatrix(k(:,jks))
end

% title
logstr = '';
if(jresponse==0)
    respstr = '';
    for j=1:length(pleGlobals.p_labels)
        respstr = [respstr ' ' pleGlobals.p_labels{j}]; %#ok<AGROW>
    end
else
    respstr = pleGlobals.p_labels{jresponse};
end
if(logplot)
    logstr = '(logarithmic)';
end
title(sprintf('response: %s %s', strrep(respstr, '_','\_'), logstr))

% save
if(exist(pleGlobals.savePath, 'dir'))
    if(isempty(jks) || length(jks) > 3)
        saveas(gcf, [pleGlobals.savePath '/relation_matrix'], 'fig')
        saveas(gcf, [pleGlobals.savePath '/relation_matrix'], 'eps')
    elseif(length(jks) == 1)
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)}], 'fig')
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)}], 'eps')
    elseif(length(jks) == 2)
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)} '_' pleGlobals.p_labels{jks(2)}], 'fig')
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)} '_' pleGlobals.p_labels{jks(2)}], 'eps')
    elseif(length(jks) == 3)
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)} '_' pleGlobals.p_labels{jks(2)} '_' pleGlobals.p_labels{jks(3)}], 'fig')
        saveas(gcf, [pleGlobals.savePath '/relation_' pleGlobals.p_labels{jks(1)} '_' pleGlobals.p_labels{jks(2)} '_' pleGlobals.p_labels{jks(3)}], 'eps')
    end
end