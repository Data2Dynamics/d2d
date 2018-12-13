% plot mcmc chains

function arPlotMCMCChains(jks, Nthinning, plot_trials)

global ar

if(~exist('Nthinning','var'))
    Nthinning = 1;
end
if(~exist('plot_trials','var'))
    plot_trials = false;
end

ps = ar.ps;
if(Nthinning>1)
    ps = ps(mod(1:size(ps,1),Nthinning)==1,:);
end
if(plot_trials)
    ps_trial = ar.ps_trial;
    if(Nthinning>1)
        ps_trial = ps_trial(mod(1:size(ps_trial,1),Nthinning)==1,:);
    end
end

h = myRaiseFigure('mcmc chains');
set(h, 'Color', [1 1 1]);

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end


jks = jks(ar.qFit(jks)==1);

[nrows, ncols] = arNtoColsAndRows(length(jks));

count = 1;
for jk=jks
    
    if(~plot_trials)
        xlimtmp2 = (max(ps(:,jk))-min(ps(:,jk)))*0.05;
        if(xlimtmp2>0)
            xlimtmp = [min(ps(:,jk))-xlimtmp2 max(ps(:,jk))+xlimtmp2];
        end
    else
        xlimtmp2 = (max(ps_trial(:,jk))-min(ps_trial(:,jk)))*0.05;
        if(xlimtmp2>0)
            xlimtmp = [min(ps_trial(:,jk))-xlimtmp2 max(ps_trial(:,jk))+xlimtmp2];
        end
    end
    
    g = subplot(nrows, ncols, count);
    arSubplotStyle(g)
    
    if(plot_trials)
        plot(ps_trial(:,jk), 'rx', 'MarkerSize', 1)
        hold on
    end
    plot(ps(:,jk), 'kx', 'MarkerSize', 1)
    hold off

    if(xlimtmp2>0)
        ylim([xlimtmp(1)-xlimtmp2*0.05 xlimtmp(2)+xlimtmp2*0.05]);
    end
    xlim([1 size(ps,1)]);
    title(arNameTrafo(ar.pLabel{jk}))
    
    
%     if(ar.xlimtmp2(count)>0)
%         ylim([ar.xlimtmp(count,1)-ar.xlimtmp2(count)*0.05 ar.xlimtmp(count,2)+ar.xlimtmp2(count)*0.05]);
%     end
%     xlim([1 size(ps,1)]);
%     title(arNameTrafo(ar.pLabel{jk}))
%     
%      ar.xlimtmp(count,:) = xlimtmp;
%      ar.xlimtmp2(count) = xlimtmp2;
     
    count = count + 1;
    
    
    
    
end



function h = myRaiseFigure(figname)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];

if(isfield(ar.ple, 'fighandel_multi') && ~isempty(ar.ple.fighandel_multi) && ...
    ar.ple.fighandel_multi ~= 0 && ...
    sum(ar.ple.fighandel_multi==openfigs)>0 && ...
    strcmp(get(ar.ple.fighandel_multi, 'Name'), figname))

    h = ar.ple.fighandel_multi;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1 0.1 0.6 0.8]);
    set(h,'Color', figcolor);
    ar.ple.fighandel_multi = h;
end


