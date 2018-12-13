% arPlotMarginalized([jks], [Nthinning])
% 
% Plots marginalized posterior of parameters in ar.ps after sampling them
% with arMC3()
% 
%   jks            parameters, which are plotted               [all]   
%   Nthinning      thinning rate for plotting                  [1]
% 
% Example:
%       arPlotMarginalized([1 2 35],10)
%
% See also arMC3, arPlotMCMCChains
function arPlotMarginalized(jks, Nthinning)

global ar

if(~exist('Nthinning','var'))
    Nthinning = 1;
end

h = myRaiseFigure('marginalized distributions');
set(h, 'Color', [1 1 1]);

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end

Nbins = 100;

jks = jks(ar.qFit(jks)==1);

[nrows, ncols] = arNtoColsAndRows(length(jks));

if(isfield(ar,'ps') && ~isempty(ar.ps))
    ps_mcmc = ar.ps;
    if(Nthinning>1)
        ps_mcmc = ps_mcmc(mod(1:size(ps_mcmc,1),Nthinning)==1,:);
    end
else
    ps_mcmc = [];
end

count = 1;
for jk=jks
    
    if(isfield(ar,'ple') && ~isempty(ar.ple) && ...
            isfield(ar.ple, 'chi2s') && ...
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
        xbins = linspace(xlimtmp(1), xlimtmp(2), Nbins);
%         ncount = histc(ps_mcmc(:,jk), xbins);
%         bar(xbins,ncount,1,'FaceColor','w');
        [ncount,xbins] = hist_norm(ps_mcmc(:,jk), xbins);
        hold on
    end
    
    % ple
    if(~isempty(ps_ple))
        % profile
        llh = transformFromLog(chi2s_ple);
        llhscale = arMatchHistWithEmpFunction(xbins, ncount, ps_ple(:,jk), llh);
        plot(ps_ple(:,jk), llh*llhscale, 'r', 'LineWidth', 1)
        hold on
        
        % thresholds
        if(ar.ple.plot_point)
            plot(xlim, transformFromLog([0 0]+min(chi2s_ple)+ar.ple.dchi2_point)*llhscale, 'r--', 'LineWidth', 1)
        end
        if(ar.ple.plot_simu)
            plot(xlim, transformFromLog([0 0]+min(chi2s_ple)+ar.ple.dchi2)*llhscale, 'r-.', 'LineWidth', 1)
        end
    end
    
    % marginalized pdf
    if(isfield(ar, 'sampling'))
        qx = ar.sampling.index==jk;
        if(sum(qx)==1)
            llhtruescale = arMatchHistWithEmpFunction(xbins, ncount, ar.sampling.ps{qx}, ar.sampling.marginalized{qx});
            plot(ar.sampling.ps{qx}, ar.sampling.marginalized{qx}*llhtruescale, 'b', 'LineWidth', 1);
        end
    end
    hold off

    xlim([xlimtmp(1)-xlimtmp2*0.05 xlimtmp(2)+xlimtmp2*0.05]);
    xlabel(['log_{10}(' arNameTrafo(ar.pLabel{jk}) ')'])
    ylabel('posterior PDF');
    
    if(count == 1)
        strleg = {};
        if(~isempty(ps_mcmc))
            strleg{end+1} = sprintf('marginalized MCMC (%i)', size(ps_mcmc,1));
        end
        if(~isempty(ps_ple))
            strleg{end+1} = 'profile posterior (PP)';
            if(ar.ple.plot_point)
                strleg{end+1} = 'PP threshold';
            end
            if(ar.ple.plot_simu)
                strleg{end+1} = 'profile likelihood threshold (simultaneous)';
            end
        end
        if(isfield(ar, 'sampling'))
            if(sum(qx)==1)
                strleg{end+1} = 'true marginalized likelihood'; %#ok<*AGROW>
            end
        end
        if(~isempty(strleg))
            legend(strleg,'Location','northoutside');
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



function b = transformFromLog(a)
b = exp(-0.5*a);

function b = transformToLog(a)
b = -2*log(a);


