% arPlotMCMCProfiles(jks, Nthinning)
% plot mcmc profiles
%
% jks       [find(ar.qFit==1)]
% Nthinning [1]

function arPlotMCMCProfiles(jks, Nthinning)

global ar

if(~exist('Nthinning','var'))
    Nthinning = 1;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

h = myRaiseFigure('mcmc profiles');
set(h, 'Color', [1 1 1]);

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end

jks = jks(ar.qFit(jks)==1);

[nrows, ncols] = arNtoColsAndRows(length(jks));

if(isfield(ar,'ps') && ~isempty(ar.ps))
    ps_mcmc = ar.ps;
    chi2s_mcmc = ar.chi2s;
    if(Nthinning>1)
        chi2s_mcmc = chi2s_mcmc(mod(1:length(chi2s_mcmc),Nthinning)==1);
        ps_mcmc = ps_mcmc(mod(1:size(ps_mcmc,1),Nthinning)==1,:);
    end
else
    ps_mcmc = [];
    chi2s_mcmc = [];
end

if(ar.config.useFitErrorMatrix==0 &&  ar.config.fiterrors == 1)
    ystr = '-2*log(L)';
    chi2s_mcmc = 2*ar.ndata*log(sqrt(2*pi)) + chi2s_mcmc;
elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
    ystr = '-2*log(L)';
    chi2s_mcmc = 2*ar.ndata_err*log(sqrt(2*pi)) + chi2s_mcmc;
else
    ystr = '\chi^2';
end

count = 1;
for jk=jks
    
    if(exist('ar.ple','var') && ~isempty(ar.ple) && ...
            ~isempty(ar.ple.chi2s) && ~isempty(ar.ple.chi2s{jk}))
        chi2s_ple = ar.ple.chi2s{jk};
        ps_ple = ar.ple.ps{jk};
    else
        chi2s_ple = [];
        ps_ple = [];
    end
    
    if(~isempty(ps_mcmc) && ~isempty(ps_ple))
        ps_all = [ps_mcmc(:,jk);ps_ple(:,jk)];
    elseif(~isempty(ps_mcmc) && isempty(ps_ple))
        ps_all = ps_mcmc(:,jk);
    elseif(isempty(ps_mcmc) && ~isempty(ps_ple))
        ps_all = ps_ple(:,jk);
    end
    
    xlimtmp2 = (max(ps_all)-min(ps_all))*0.05;
    if(xlimtmp2>0)
        xlimtmp = [min(ps_all)-xlimtmp2 max(ps_all)+xlimtmp2];
    end
    
    g = subplot(nrows, ncols, count);
    arSubplotStyle(g);
    
    % plot MCMC
    if(~isempty(ps_mcmc))
        plot(ps_mcmc(:,jk), chi2s_mcmc, 'bx')
        hold on
    end
    
    % ple
    if(~isempty(ps_ple))
        % profile
        plot(ps_ple(:,jk), chi2s_ple, 'k', 'LineWidth', 1)
        hold on
        
        % thresholds
        if(ar.ple.plot_point)
            plot(xlim, [0 0]+min(chi2s_ple)+ar.ple.dchi2_point, 'r--', 'LineWidth', 1)
        end
        if(ar.ple.plot_simu)
            plot(xlim, [0 0]+min(chi2s_ple)+ar.ple.dchi2, 'r-.', 'LineWidth', 1)
        end
    end
    hold off

    xlim([xlimtmp(1)-xlimtmp2*0.05 xlimtmp(2)+xlimtmp2*0.05]);
    xlabel(['log_{10}(' arNameTrafo(ar.ple.p_labels{jk}) ')'])
    ylabel(ystr);
    
    if(~isempty(ps_ple))
        if(ar.ple.plot_simu)
            ylim([min(chi2s_ple) min(chi2s_ple)+1.2*ar.ple.dchi2])
        else
            ylim([min(chi2s_ple) min(chi2s_ple)+1.2*ar.ple.dchi2_point])
        end
    end
    
    if(count == 1)
        strleg = {};
        if(~isempty(ps_mcmc))
            strleg{end+1} = sprintf('samples (%i)', length(chi2s_mcmc));
        end
        if(~isempty(ps_ple))
            strleg{end+1} = 'profile likelihood';
            if(ar.ple.plot_point)
                strleg{end+1} = 'profile likelihood threshold (point-wise)';
            end
            if(ar.ple.plot_simu)
                strleg{end+1} = 'profile likelihood threshold (simultaneous)';
            end
        end
        if(~isempty(strleg))
            legend(strleg);
        end
    end
    
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


