% Plot compact profile Likelihoods
%
% plePlotMulti(jks, savetofile, ncols, nrows, show_hit_bound, plot_hit_bound, plot_thresholds)
%
% indices           plot only parameters jks                        [all]
% savetofile        save plot                                       [false]
% ncols
% nrows
% show_hit_bound    index show hitting boundary of parameters       [all]
% plot_hit_bound    highlight hitting boundary of parameters        [true]
% plot_thresholds   plot confidence threshold                       [true]

function plePlotMulti(jks, savetofile, ncols, nrows, show_hit_bound, plot_hit_bound, plot_thresholds)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(~isfield(ar.ple,'ps') || isempty(ar.ple.ps))
    return
end
if(~exist('jks','var') || isempty(jks))
    jks = 1:length(ar.p);
end
if(~exist('savetofile','var') || isempty(savetofile))
    savetofile = false;
end

% Map PLE to current ar struct
if(ischar(jks))
    jks = find(ismember(ar.ple.p_labels, jks)); %R2013a compatible
else
    jks = find(ismember(ar.ple.p_labels, ar.pLabel(jks))); %R2013a compatible
end
pleToP = ismember(ar.pLabel, ar.ple.p_labels);
pToPLE = ismember(ar.ple.p_labels, ar.pLabel);
lb(pToPLE) = ar.lb(pleToP);
ub(pToPLE) = ar.ub(pleToP);
qFit(pToPLE) = ar.qFit(pleToP);
qLog10(pToPLE) = ar.qLog10(pleToP);
jks(jks>numel(ar.ple.ps)) = [];

sumples = 0;
for j=jks
    if(~isempty(ar.ple.ps{j}))
        sumples = sumples + 1;
    end
end

if(~exist('ncols', 'var') || isempty(ncols))
    ncols = ceil(sumples^(0.4))+1;
    nrows = ceil(sumples/ncols);
end
if(~exist('nrows', 'var') || isempty(nrows))
    nrows = ceil(sumples/ncols);
end
if(~exist('show_hit_bound','var') || isempty(show_hit_bound))
    show_hit_bound = 1:length(ar.ple.ps);
end
if(~exist('plot_hit_bound','var') || isempty(plot_hit_bound))
    plot_hit_bound = true;
end
if(~exist('plot_thresholds','var') || isempty(plot_thresholds))
    plot_thresholds = true;
end

if(ar.ple.plot_point && ~ar.ple.plot_simu)
    strtitle = sprintf('profile log-likelihood (point-wise)');
elseif(~ar.ple.plot_point && ar.ple.plot_simu)
    strtitle = sprintf('profile log-likelihood (simultaneous)');
else
    strtitle = sprintf('profile log-likelihood');
end

h = myRaiseFigure(strtitle);
set(h, 'Color', [1 1 1]);

if ( isfield( ar.ple, 'firstID' ) )
    count = ar.ple.firstID;
else
    count = 1;
end

minchi2 = Inf;
for jk=jks
    if(~isempty(ar.ple.ps{jk}))
        minchi2 = min([minchi2 min(ar.ple.chi2s{jk})]);
    end
end

for jk=jks
    if(~isempty(ar.ple.ps{jk}))
        g = subplot(nrows,ncols,count);
        arSubplotStyle(g);
        
        ps = ar.ple.ps{jk};
        chi2s = ar.ple.chi2s{jk};
        
        qCloseToUB = ps > ones(length(chi2s),1) * (ub - ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * qFit==1;
        qCloseToLB = ps < ones(length(chi2s),1) * (lb + ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * qFit==1;
     
        qhitbound = false(size(ps));
        qhitbound(:,qFit==1) = ar.ple.gradient{jk}(:,qFit==1) > ar.ple.grad_thres & qCloseToLB(:,qFit==1) | ...
            ar.ple.gradient{jk}(:,qFit==1) < -ar.ple.grad_thres & qCloseToUB(:,qFit==1);
        
        % profile
        qplot = true(size(ps,1),1);
        if(sum(show_hit_bound==jk)>0)
            plot(ps(:,jk), chi2s, 'k-', 'LineWidth', 1)
        else
            qplot = sum(qhitbound,2)==0;
            plot(ps(qplot,jk), chi2s(qplot), 'k-', 'LineWidth', 1)
        end
        hold on
        if(isfield(ar.ple,'chi2spriors'))
            mod_const = min(chi2s(qplot) - ar.ple.chi2spriors{jk}(qplot));
            plot(ps(qplot,jk), ar.ple.chi2spriors{jk}(qplot) + mod_const, 'b--', 'LineWidth', 1)
        end
        if(isfield(ar.ple,'chi2spriorsAll'))
            mod_const = min(chi2s(qplot) - ar.ple.chi2spriorsAll{jk}(qplot));
            plot(ps(qplot,jk), ar.ple.chi2spriorsAll{jk}(qplot) + mod_const, 'c--', 'LineWidth', 1)
        end
        
        % boundary values
        if(sum(show_hit_bound==jk)>0 && plot_hit_bound)
            plot(ps(sum(qhitbound,2)>0,jk), chi2s(sum(qhitbound,2)>0), 'ko', 'LineWidth', 1)
        end
        
        % limits
        dchi2 = ar.ple.dchi2_point;
        if(ar.ple.plot_simu)
            dchi2 = ar.ple.dchi2;
        end

        ylimmax = ar.ple.merit+dchi2*1.3;
        if(plot_thresholds)
            ylim([minchi2-dchi2*0.1 ylimmax]);
        end
        
        qbelowchi2 = chi2s < ylimmax;
        xlimtmp2 = (max(ps(qbelowchi2,jk))-min(ps(qbelowchi2,jk)))*0.05;
        if(xlimtmp2>0)
            if(show_hit_bound)
                xlimtmp = [min(ps(qbelowchi2,jk))-xlimtmp2 max(ps(qbelowchi2,jk))+xlimtmp2];
            else
                xlimtmp = [min(ps(sum(qhitbound,2)==0 & qbelowchi2,jk))- ...
                    xlimtmp2 max(ps(sum(qhitbound,2)==0 & qbelowchi2,jk))+xlimtmp2];
            end
            xlim(xlimtmp);
        end
        
        % thresholds
        if(plot_thresholds)
            if(ar.ple.plot_point)
                plot(xlim, [0 0]+minchi2+chi2inv(1-ar.ple.alpha_level, 1), 'r--')
            end
            if(ar.ple.plot_simu)
                plot(xlim, [0 0]+minchi2+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), 'r--')
            end
            
            if(count == 1)
                if(ar.ple.plot_point && ~ar.ple.plot_simu)
                    text(mean(xlim), minchi2+chi2inv(1-ar.ple.alpha_level, 1), sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                elseif(~ar.ple.plot_point && ar.ple.plot_simu)
                    text(mean(xlim), minchi2+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                else
                    text(mean(xlim), minchi2+chi2inv(1-ar.ple.alpha_level, 1), sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                    text(mean(xlim), minchi2+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
                end
            end
        end
        
        % hessian CI
%         sampletmp = linspace(ar.ple.conf_lb(1,jk), ar.ple.conf_ub(1,jk), 100);
%         assym_chi2 = min(chi2s) + (1/ar.ple.pstd(jk)^2)*(ar.ple.p(jk)-sampletmp).^2;
%         if(isreal(sampletmp))
%             plot(sampletmp, assym_chi2, '-', 'Color', [.5 .5 .5], 'LineWidth', 1)
%         end
        
        % optimum
        plot(ar.ple.p(jk), ar.ple.merit, '*', 'Color', [.5 .5 .5], 'LineWidth', 1)
        hold off

        if ( qLog10(jk) )
            xlabel(['log_{10}(' arNameTrafo(ar.ple.p_labels{jk}) ')']);
        else
            xlabel(arNameTrafo(ar.ple.p_labels{jk}));
        end
        title(sprintf('parameter #%i', find(strcmp(ar.pLabel, ar.ple.p_labels{jk}))));
        
        if(mod(count-1,ncols)==0)
            ylabel(ar.ple.ylabel)
        else
            ylabel('')
            if(plot_thresholds)
                set(gca, 'YTickLabel', {})
            end
        end
        
        count = count + 1;
    end
end

% save

if savetofile
    if ~exist(ar.ple.savePath, 'dir')
        mkdir(ar.ple.savePath)
    end
    
    name = 'multi_plot';
    ar.ple.figPathMulti{jk} = [ar.ple.savePath '/' name];
    set(h,'Renderer','painters')
    savePath = [ar.ple.savePath '/' name];
    saveas(h, savePath, 'fig');
    
    
    if(ispc)
        print('-dpdf', savePath);
    elseif(isunix)
        set(h,'Units','in')
        myaxes = findobj(h,'Type','axes');
        mypos = get(myaxes,'Position');
        if(iscell(mypos))
            mypos = cell2mat(mypos);
        end
        mycols = unique(round(mypos(:,1)*1e6)/1e6);
        myrows = unique(round(mypos(:,2)*1e6)/1e6);
        nCols = length(mycols);
        nRows = length(myrows);

        myfigpos = get(h,'Position');
        set(h,'Position',[myfigpos(1:2) 4*nCols 4*nRows]);

        for ia = 1:length(myaxes)
            [tmp, indCol] = min(abs(mypos(ia,1) - mycols));
            [tmp, indRow] = min(abs(mypos(ia,2) - myrows));

            set(myaxes(ia),'Position',[(indCol-1)/nCols+(1/nCols)*0.15 (indRow-1)/nRows+(1/nRows)*0.1 (1/nCols)*0.7 (1/nRows)*0.8])
        end
        axoptions={'x tick label style={/pgf/number format/fixed}','y tick label style={/pgf/number format/fixed}'};    
        matlab2tikz([savePath '_Report.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'width','0.9\textwidth','automaticLabels',true,'extraAxisOptions',axoptions)
        matlab2tikz([savePath '.tex'],'figurehandle',h,'showInfo', false, 'showWarnings',false, 'standalone', true,'automaticLabels',true,'extraAxisOptions',axoptions)

        library_path = getenv('LD_LIBRARY_PATH');
        setenv('LD_LIBRARY_PATH', '');
        
        workingdir = cd;
        cd(ar.ple.savePath);

        if(ismac)
            eval(['!/usr/texbin/pdflatex -interaction=nonstopmode ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
        else
            eval(['!pdflatex -interaction=nonstopmode ' [arPathConvert(name) '.tex'] ' > log_pdflatex.txt']);
        end

        cd(workingdir);
        setenv('LD_LIBRARY_PATH', library_path);

        for ia = 1:length(myaxes)
            set(myaxes(ia),'Position',mypos(ia,:))
        end
        set(h,'Position',myfigpos);
    end
    
    if(which('plot2svg'))
%         plot2svg([savePath '.svg'],h)
    end
    
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


