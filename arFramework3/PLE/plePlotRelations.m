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

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
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
    
    for j=find(ar.qFit)
        if (j<=size(ar.ple.ps,2))
            if(~isempty(ar.ple.ps{j}))
                k = [k;ar.ple.ps{j}]; %#ok<AGROW>
                chi2 = [chi2 ar.ple.chi2s{j}]; %#ok<AGROW>
            end
        end
    end
else
    k = ar.ple.ps{jresponse};
    chi2 = ar.ple.chi2s{jresponse};
end

if(logplot)
    k = log10(k);
end

nfarben = 50;

% farben = winter(nfarben);
% farben(end,:) = [1 0 0];

farben = gray(nfarben);
farben(end,:) = [1 1 1];

farbenindex = (chi2-min(chi2))./ar.ple.dchi2;
farbenindex = 1+floor(farbenindex*nfarben);
farbenindex(farbenindex>nfarben) = nfarben;

figure(length(ar.ple.p)+2)

if(isempty(jks))
    plotmatrix(k(:,ar.qFit))
elseif(length(jks) == 1)
    hist(k(:,jks), ...
        sum(~isnan(k(:,jks)))/5)
    xlabel(ar.ple.p_labels{jks(1)})
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
    xlabel(ar.ple.p_labels{jks(1)})
    ylabel(ar.ple.p_labels{jks(2)})
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
    xlabel(ar.ple.p_labels{jks(1)})
    ylabel(ar.ple.p_labels{jks(2)})
    zlabel(ar.ple.p_labels{jks(3)})
    grid on
elseif(length(jks) > 3)
    plotmatrix(k(:,jks))
end

% title
logstr = '';
if(jresponse==0)
    respstr = '';
    for j=1:length(ar.ple.p_labels)
        respstr = [respstr ' ' ar.ple.p_labels{j}]; %#ok<AGROW>
    end
else
    respstr = ar.ple.p_labels{jresponse};
end
if(logplot)
    logstr = '(logarithmic)';
end
title(sprintf('response: %s %s', strrep(respstr, '_','\_'), logstr))

% save
if(exist(ar.ple.savePath, 'dir'))
    if(isempty(jks) || length(jks) > 3)
        saveas(gcf, [ar.ple.savePath '/relation_matrix'], 'fig')
        saveas(gcf, [ar.ple.savePath '/relation_matrix'], 'eps')
    elseif(length(jks) == 1)
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)}], 'fig')
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)}], 'eps')
    elseif(length(jks) == 2)
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)} '_' ar.ple.p_labels{jks(2)}], 'fig')
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)} '_' ar.ple.p_labels{jks(2)}], 'eps')
    elseif(length(jks) == 3)
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)} '_' ar.ple.p_labels{jks(2)} '_' ar.ple.p_labels{jks(3)}], 'fig')
        saveas(gcf, [ar.ple.savePath '/relation_' ar.ple.p_labels{jks(1)} '_' ar.ple.p_labels{jks(2)} '_' ar.ple.p_labels{jks(3)}], 'eps')
    end
end