function [] = plotMulti_PESTO(ar)
    %% SORT RESULTS and COMPUTE OBJECTIVES
    xx = sort(ar.chi2s);
    xx_logL = -ar.ndata*log(sqrt(2*pi)) - xx./2;
    xx_logL = xx_logL';

    %% CLUSTERING
    n_starts = length(xx_logL);
    if (n_starts > 1)
        clust = cluster(linkage(pdist(xx_logL)),'cutoff',0.1,'criterion','distance');
    else
        clust = 1;
    end
    uclust = unique(clust);

    for iclust = 1:length(uclust)
        sizecluster(iclust) = sum(clust == uclust(iclust));
    end

    %% ASSIGN COLORS
    Col = colormap(gray(n_starts+ceil(n_starts/3)));
    Col = Col.^(1/3);
    Col(1,:) = [1,0,0];

    % sort clusters
    for iclust = 1:length(uclust)
        Jclust(iclust) = max(xx_logL(find(clust == uclust(iclust))));
    end
    Jclust(isnan(Jclust)) = -Inf;
    [~,idx] = sort(Jclust,'descend');
    uclust = uclust(idx);
    sizecluster = sizecluster(idx);

    if(sizecluster(1)>1)
        ColClust = [1,0,0;flipud(parula(max(sum(sizecluster>1)-1,0)))];
    else
        ColClust = flipud(parula(sum(sizecluster>1)));
    end

    for iclust = 1:length(uclust)
        if(sizecluster(iclust)>1)
        Col(clust == uclust(iclust),:) = repmat(ColClust(sum(sizecluster(1:iclust)>1),:),[sizecluster(iclust),1]);
        end
    end

    figure;
    
    %% PLOT OBJECTIVES
    subplot(2,2,1);
    n_finished_starts = 0;
    for j = 1 : n_starts
        if ~isnan(xx_logL(j))
            n_finished_starts = j;
        else
            break;
        end
    end

    plot(1:n_finished_starts,xx_logL(1:n_finished_starts),'-','color',0.9*[1,1,1],'linewidth',2);
    hold on;
    for j = n_finished_starts:-1:1
        plot(j,xx_logL(j),'o','color',Col(j,:),'linewidth',2);
        hold on;
    end
    hold off;
    xlim([1-0.2,n_starts+0.2]);
    xlabel('start');
    ylabel('log-likelihood');

    title('all estimates');


    %% PLOT TOP TEN OBJECTIVES
    subplot(2,2,3);
    plot(1:min(n_starts,10),xx_logL(1:min(n_starts,10)),'-','color',0.9*[1,1,1],'linewidth',2); hold on;
    for j = min(n_starts,10):-1:1
        plot(j,xx_logL(j),'o','color',Col(j,:),'linewidth',2); hold on;
    end
    hold off;
    xlim([1-0.2,min(n_starts,10)+0.2]);
    if(any(~isnan(xx_logL(1:min(n_starts,10)))))
        ylim([min(xx_logL(1),min(xx_logL(1:min(n_starts,10)))-1),xx_logL(1)+1]);
    end
    xlabel('start');
    ylabel('log-likelihood');

    title('top 10 estimates');
end