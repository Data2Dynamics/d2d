% plot scan of likelihood

function arPlotScan(jks)

global ar

if(~exist('jks','var') || isempty(jks))
    jks = find(~cellfun(@isempty, ar.scan.ps));
end

sumples = 0;
for j=jks
    if(~isempty(ar.scan.ps{j}))
        sumples = sumples + 1;
    end
end

[nrows, ncols] = arNtoColsAndRows(sumples);


h = arRaiseFigure(ar, 'fighandel_multi', 'likelihood scan');
clf
set(h, 'Color', [1 1 1]);

count = 1;

minchi2 = Inf;
for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        minchi2 = min([minchi2 min(ar.scan.chi2s{jk})]);
    end
end
chi2curr = ar.chi2fit;
if(ar.config.fiterrors == 1)
    chi2curr = 2*ar.ndata*log(sqrt(2*pi)) + chi2curr;
end

for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        g = subplot(nrows,ncols,count);
        
        ps = ar.scan.ps{jk};
        chi2s = ar.scan.chi2s{jk};
        constrs = ar.scan.constrs{jk};
        
        if(ar.config.fiterrors == 1)
            chi2s = 2*ar.ndata*log(sqrt(2*pi)) + chi2s;
        end
        
        if(ar.ndata>0)
            plot(ps, chi2s, 'k-', 'LineWidth', 1)
            hold on
        end
        if(ar.nconstr>0)
            plot(ps, constrs, 'r--', 'LineWidth', 1)
        end
        hold off
        
        ax1 = g;
        
        % optimum
        line(ar.p(jk), chi2curr, 'Marker', '*', 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerSize', 8, ...
            'Parent', ax1)
        

        xlabel(ax1, ['log_{10}(' arNameTrafo(ar.pLabel{jk}) ')'])
        arSpacedAxisLimits(g)
        
        if(mod(count-1,ncols)==0)
            if(ar.config.fiterrors == 1)
                ylabel(ax1, '-2*log(L)');
            else
                ylabel(ax1, '\chi^2');
            end
        end
        
        count = count + 1;
    end
end
