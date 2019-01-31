% plot one-dimensional cut of likelihood.
% All other parameter are fixed
% 
% The function is intended to be run after arScan.m

function arPlotScan(jks, savetofile)

global ar

if ~isfield(ar,'scan')
    warning('ar.scan not available. Run arScan first1')
    return
end

if(~exist('jks','var') || isempty(jks))
    jks = find(~cellfun(@isempty, ar.scan.ps));
end
if(~exist('savetofile','var'))
    savetofile = false;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

sumples = 0;
for j=jks
    if(~isempty(ar.scan.ps{j}))
        sumples = sumples + 1;
    end
end

[nrows, ncols] = arNtoColsAndRows(sumples);


h = myRaiseFigure('likelihood scan');
clf
set(h, 'Color', [1 1 1]);

count = 1;

minchi2 = Inf;
for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        minchi2 = min([minchi2 min(ar.scan.chi2s{jk})]);
    end
end
chi2curr = arGetMerit;
% if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
%     chi2curr = 2*ar.ndata*log(sqrt(2*pi)) + chi2curr;
% elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
%     chi2curr = 2*ar.ndata_err*log(sqrt(2*pi)) + chi2curr;
% end

for jk=jks
    if(~isempty(ar.scan.ps{jk}))
        g = subplot(nrows,ncols,count);
        
        ps = ar.scan.ps{jk};
        chi2s = ar.scan.chi2s{jk};
        constrs = ar.scan.constrs{jk};
        
%         if(ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1)
%             chi2s = 2*ar.ndata*log(sqrt(2*pi)) + chi2s;
%         elseif(ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0)
%             chi2s = 2*ar.ndata_err*log(sqrt(2*pi)) + chi2s;
%         end
        
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
            if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors == 1) || ...
                    (ar.config.useFitErrorMatrix==1 && sum(sum(ar.config.fiterrors_matrix==1))>0) )
                ylabel(ax1, '-2*log(L)');
            else
                ylabel(ax1, '\chi^2');
            end
        end
        
        count = count + 1;
    end
end

% save
if(savetofile && exist(ar.ple.savePath, 'dir'))
    ar.ple.figPathMulti{jk} = [ar.ple.savePath '/multi_plot'];
    save([ar.ple.savePath '/results.mat'], 'ar.ple');
    saveas(gcf, [ar.ple.savePath '/multi_plot'], 'fig')
    print('-depsc2', [ar.ple.savePath '/multi_plot']);
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
