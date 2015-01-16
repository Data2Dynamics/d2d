% plot mcmc chains

function arPlotMCMCParVsChi(jks, Nthinning, popt, chi2opt)

global ar

if(~exist('Nthinning','var'))
    Nthinning = 1;
end
ps = ar.ps;
chi2s = ar.chi2s;
if(Nthinning>1)
    ps = ps(mod(1:size(ps,1),Nthinning)==1,:);
    chi2s = chi2s(mod(1:length(chi2s),Nthinning)==1,:);
end

h = arRaiseFigure(pleGlobals, 'fighandel_multi', 'mcmc chains');
set(h, 'Color', [1 1 1]);

if(~exist('jks','var') || isempty(jks))
    jks = find(ar.qFit==1);
end


jks = jks(ar.qFit(jks)==1);

[nrows, ncols] = arNtoColsAndRows(length(jks));

count = 1;
for jk=jks
    
    xlimtmp2 = (max(ps(:,jk))-min(ps(:,jk)))*0.05;
    if(xlimtmp2>0)
        xlimtmp = [min(ps(:,jk))-xlimtmp2 max(ps(:,jk))+xlimtmp2];
    end
    
    g = subplot(nrows, ncols, count);
    arSubplotStyle(g)

    plot(ps(:,jk), chi2s, 'kx', 'MarkerSize', 1)
    if(nargin>2)
        hold on
        % colors = jet(length(popt));
        for j=1:length(popt)
            % plot(popt{j}(jk), chi2opt(j), '*', 'Color', colors(j,:));
            plot(popt{j}(jk), chi2opt(j), '*r');
        end
        hold off
    end
    
    xlim([xlimtmp(1)-xlimtmp2*0.05 xlimtmp(2)+xlimtmp2*0.05]);
    
    if(nargin>2)
        ylim([min([chi2s(:); chi2opt(:)]) quantile(chi2s, 0.99)]);
    else
        ylim([min(chi2s(:)) quantile(chi2s, 0.99)]);
    end
    title(arNameTrafo(ar.pLabel{jk}))
    
    count = count + 1;
end



