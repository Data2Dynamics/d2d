% Plot profile likelohood
%
% h = plePlot(indices, savetofile, subs_para, para_index, show_hit_bound, plot_hit_bound)
%
% indices           plot only parameters jks                        [all]
% savetofile        save plot                                       [false]
% subs_para         show parameter values more nice                 [true]
% para_index        plot only relations to parameters para_index    [all]
% show_hit_bound    show hitting boundary of parameters             [true]
% plot_hit_bound    plot hitting boundary of parameters             [true]
% 
%   h   handle to figure object
function varargout = plePlot(indices, savetofile, subs_para, para_index, show_hit_bound, plot_hit_bound)

global ar

if(isempty(ar.ple))
    error('perform ple before usage');
end
if(isempty(ar.ple.ps))
    return
end

if(~exist('indices','var') || isempty(indices))
    indices = 1:length(ar.ple.ps);
elseif ischar(indices)
    indices = strmatch(indices,ar.ple.p_labels,'exact');
elseif ischar(indices)
    [~,indices] = intersect(ar.ple.p_labels,indices);
end
if(~exist('savetofile','var'))
    savetofile = false;
end
if(~exist('subs_para','var'))
    subs_para = true;
end
if(~exist('show_hit_bound','var'))
    show_hit_bound = true;
end
if(~exist('plot_hit_bound','var'))
    plot_hit_bound = true;
end

if ( savetofile )
    if(~exist(ar.ple.savePath, 'dir'))
        mkdir('.', ar.ple.savePath)
    end
end

farben = lines(length(ar.ple.p_labels));
zeichen = {'-', '--', '-.'};
ps_label_count = 5;
ps_var_min_label = 1e-2;

if(~isfield(ar.ple, 'fighandel'))
    ar.ple.fighandel = zeros(size(ar.ple.p));
end

count = 0;
for jj=1:length(indices)
    jk = indices(jj);
    if(~isempty(ar.ple.ps{jk}))
        h = myRaiseFigure(jk, sprintf('PLE#%i %s', jk, ar.ple.p_labels{jk}), count);
        set(h, 'Color', [1 1 1]);
        
        g = subplot(5,1,[1 2 3]);
        arSubplotStyle(g);
        
        ps = ar.ple.ps{jk};
        chi2s = ar.ple.chi2s{jk};
        
        qCloseToUB = ps > ones(length(chi2s),1) * (ar.ub - ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * ar.qFit==1;
        qCloseToLB = ps < ones(length(chi2s),1) * (ar.lb + ar.ple.dist_thres) & ...
            ones(length(chi2s),1) * ar.qFit==1;

        qhitbound = false(size(ps));
        qhitbound(:,ar.qFit==1) = ar.ple.gradient{jk}(:,ar.qFit==1) > ar.ple.grad_thres & qCloseToLB(:,ar.qFit==1) | ...
            ar.ple.gradient{jk}(:,ar.qFit==1) < -ar.ple.grad_thres & qCloseToUB(:,ar.qFit==1);
        
        % profile
        if(show_hit_bound)
            plot(ps(:,jk), chi2s, 'k-', 'LineWidth', 1)
            hold on
            if(isfield(ar.ple, 'chi2sviolations'))
                plot(ps(:,jk), chi2s-ar.ple.chi2sviolations{jk}, 'b', 'LineWidth', 1)
            end
            if(isfield(ar.ple,'chi2spriors'))
                mod_const = min(chi2s - ar.ple.chi2spriors{jk});
                plot(ps(:,jk), ar.ple.chi2spriors{jk} + mod_const, 'b--', 'LineWidth', 1)
            end
        else
            qplot = sum(qhitbound,2)==0;
            plot(ps(qplot,jk), chi2s(qplot), 'k-', 'LineWidth', 1)
            hold on
%             if(isfield(ar.ple, 'chi2sviolations'))
%                 plot(ps(qplot,jk), chi2s(qplot)-ar.ple.chi2sviolations{jk}(qplot), 'b', 'LineWidth', 1)
%             end
            if(isfield(ar.ple,'chi2spriors'))
                mod_const = min(chi2s(qplot) - ar.ple.chi2spriors{jk}(qplot));
                plot(ps(qplot,jk), ar.ple.chi2spriors{jk}(qplot) + mod_const, 'b--', 'LineWidth', 1)
            end
        end
        
        % boundary values
        if(show_hit_bound && plot_hit_bound)
            plot(ps(sum(qhitbound,2)>0,jk), chi2s(sum(qhitbound,2)>0), 'ko', 'LineWidth', 1)
        end
        
        % limits
        dchi2 = ar.ple.dchi2_point;
        if(ar.ple.plot_simu)
            dchi2 = ar.ple.dchi2;
        end
        
        ylimmax = min(chi2s)+dchi2*1.3;
%         if(isfield(ar.ple, 'chi2sviolations'))
%             ylim([min(chi2s-ar.ple.chi2sviolations{jk})-dchi2*0.1 ylimmax]);
%         else
            ylim([min(chi2s)-dchi2*0.1 ylimmax]);
%         end
        
        qbelowchi2 = chi2s < ylimmax;
        xlimtmp2 = (max(ps(qbelowchi2,jk))-min(ps(qbelowchi2,jk)))*0.05;
        if(xlimtmp2>0)
            if(show_hit_bound)
                xlimtmp = [min(ps(qbelowchi2,jk))-xlimtmp2 max(ps(qbelowchi2,jk))+xlimtmp2];
            else
                xlimtmp = [min(ps(sum(qhitbound,2)==0 & qbelowchi2',jk))- ...
                    xlimtmp2 max(ps(sum(qhitbound,2)==0 & qbelowchi2',jk))+xlimtmp2];
            end
            xlim(xlimtmp);
        end
        
        % thresholds
        if(ar.ple.plot_point)
            plot(xlim, [0 0]+min(chi2s)+chi2inv(1-ar.ple.alpha_level, 1), 'r--')
        end
        if(ar.ple.plot_simu)
            plot(xlim, [0 0]+min(chi2s)+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), 'r--')
        end
        
        if(ar.ple.plot_point && ~ar.ple.plot_simu)
            text(mean(xlim), min(chi2s)+chi2inv(1-ar.ple.alpha_level, 1), sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        elseif(~ar.ple.plot_point && ar.ple.plot_simu)
            text(mean(xlim), min(chi2s)+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        else
            text(mean(xlim), min(chi2s)+chi2inv(1-ar.ple.alpha_level, 1), sprintf('%2i%% (point-wise)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
            text(mean(xlim), min(chi2s)+chi2inv(1-ar.ple.alpha_level, ar.ple.dof), sprintf('%2i%% (simultaneous)', (1-ar.ple.alpha_level)*100), 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
            
        % ple CIs
%         if(ar.ple.plot_simu)
%             plot([0 0]+ar.ple.conf_lb(2,jk), ylim, '-', 'Color', [1 0 0])
%             plot([0 0]+ar.ple.conf_ub(2,jk), ylim, '-', 'Color', [1 0 0])
%         end
%         if(ar.ple.plot_point)
%             plot([0 0]+ar.ple.conf_lb_point(2,jk), ylim, '--', 'Color', [1 0 0])
%             plot([0 0]+ar.ple.conf_ub_point(2,jk), ylim, '--', 'Color', [1 0 0])
%         end

        % plot parameters bounds
%         plot([0 0]+ar.lb(jk), ylim, 'k--')
%         plot([0 0]+ar.ub(jk), ylim, 'k--')

        % plot true value
%         if(isfield(ar.ple, 'p_true'))
%             plot([0 0]+ar.ple.p_true(jk), ylim, 'g--')
%         end
        
        % hessian CI
%         xlimtmp = xlim;
%         sampletmp = linspace(xlimtmp(1), xlimtmp(2), 100);
%         assym_chi2 = min(chi2s) + (1/ar.ple.pstd(jk)^2)*(ar.ple.p(jk)-sampletmp).^2;
%         if(isreal(sampletmp))
%             plot(sampletmp, assym_chi2, '-.', 'Color', [.5 .5 .5])
%         end

        % optimum
        plot(ar.ple.p(jk), ar.ple.merit, '*', 'Color', [.5 .5 .5], 'LineWidth', 1, 'MarkerSize', 8)
        hold off
        
        ylabel(ar.ple.ylabel)
        
        % parameter changes
        
        if(~exist('para_index','var') || isempty(para_index))
            notjk = 1:length(ar.ple.p);
            notjk = notjk~=jk;
            notjk = notjk & ar.qFit==1;
        else
            notjk = zeros(size(ar.ple.p));
            notjk(para_index) = 1;
        end
        
        legendstmp = cell(1,sum(notjk));
        legendstmplines = zeros(1,sum(notjk));
        ps_var = zeros(1,sum(notjk));
        
        g = subplot(5,1,[4 5]);
        arSubplotStyle(g);
        
        if(subs_para==2)
            ps(:,notjk) = bsxfun(@minus,ps(:,notjk),ar.lb(notjk));
            ps(:,notjk) = bsxfun(@rdivide,ps(:,notjk),ar.ub(notjk)-ar.lb(notjk));
        end
        
        ccount = 1;
        for j=find(notjk)
            qisnonan = ~isnan(ps(:,j));
            ps_var(ccount) = max(ps(qisnonan,j)) - min(ps(qisnonan,j));
            zeichenindex = mod(floor((j-1)/7)+1, 3)+1;

            if(subs_para==1)
                medianp = ar.ple.p(j);
%                 medianp = median(ps(qisnonan,j));
            else
                medianp = 0;
            end
            
            if(show_hit_bound)
                line_s = plot(ps(:,jk), ps(:,j)-medianp, [zeichen{zeichenindex}], 'Color', farben(j,:));
            else
                line_s = plot(ps(sum(qhitbound,2)==0,jk), ps(sum(qhitbound,2)==0,j)-medianp, [zeichen{zeichenindex}], 'Color', farben(j,:));
            end
            hold on
            if(show_hit_bound && plot_hit_bound)
                plot(ps(qhitbound(:,j),jk), ps(qhitbound(:,j),j)-medianp, 'o', 'Color', farben(j,:));
            end
            
            legendstmp{ccount} = arNameTrafo(ar.ple.p_labels{j}); %#ok<*AGROW>
            legendstmplines(ccount) = line_s;
            ccount = ccount + 1;
        end 
        hold off
        if(~isempty(legendstmplines))
            if(ps_label_count<length(legendstmp))
                [ps_varsorted, i_largest_std] = sort(ps_var, 2, 'descend');
                i_largest_std = i_largest_std(1:ps_label_count);
                qi_largest_std = ps_varsorted(1:ps_label_count) > ps_var_min_label;
                i_largest_std = i_largest_std(qi_largest_std);
                if(~isempty(i_largest_std))
                    legend(legendstmplines(i_largest_std), legendstmp(i_largest_std), 'Location', 'Best')
                end
            else
                legend(legendstmplines, legendstmp, 'Location', 'Best')
            end
        end
        
        % If the interactivity system is enabled, register the callbacks
        % and provide arInteractivity with the required data.
        if ( exist( 'arInteractivity', 'file' ) )
            if ( arInteractivity )
                if (~exist('i_largest_std', 'var')|~i_largest_std) %#ok
                    i_largest_std = 1 : length( legendstmp );
                end
                lineLegends.legends = legendstmp;
                lineLegends.handles = legendstmplines;
                lineLegends.currentLegend.handles = legendstmplines(i_largest_std);
                lineLegends.currentLegend.legends = legendstmp(i_largest_std);
                lineLegends.currentLegend.ID = 1;
                arInteractivity( 'ple', lineLegends );
            end
        end
        
        if(xlimtmp2>0)
            xlim(xlimtmp);
        end
        ylabel({'change of';'other parameters'})
        if(ar.qLog10(jk))
            xlabel(['log_{10}(' arNameTrafo(ar.ple.p_labels{jk}) ')'])
        else
            xlabel(arNameTrafo(ar.ple.p_labels{jk}))
        end
        grid(g,'on');
%         axis(g,'square');
        
        % save
        if(savetofile && exist(ar.ple.savePath, 'dir'))
            ar.ple.figPath{jk} = [ar.ple.savePath '/' ar.ple.p_labels{jk}];
            saveas(gcf, [ar.ple.savePath '/' ar.ple.p_labels{jk}], 'fig')
            print('-depsc2', [ar.ple.savePath '/' ar.ple.p_labels{jk}]);
            if(ispc)
                print('-dpdf', [ar.ple.savePath '/' ar.ple.p_labels{jk}]);
            elseif(ismac)
                system(['ps2pdf  -dEPSCrop ' [ar.ple.savePath '/' ar.ple.p_labels{jk}] '.eps '  [ar.ple.savePath '/' ar.ple.p_labels{jk}] '.pdf']);
            else
                system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' [ar.ple.savePath '/' ar.ple.p_labels{jk}] '.eps '  [ar.ple.savePath '/' ar.ple.p_labels{jk}] '.pdf']);
            end
        end
        
        count = count + 1;
    end
end

if nargout>0
    varargout{1} = h;
end



function h = myRaiseFigure(jk, figname, figcount)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.01;

if(isfield(ar.ple, 'fighandel') && ~isempty(ar.ple.fighandel) && ...
    ar.ple.fighandel(jk) ~= 0 && ...
    sum(ar.ple.fighandel(jk)==openfigs)>0 && ...
    strcmp(get(ar.ple.fighandel(jk), 'Name'), figname))

    h = ar.ple.fighandel(jk);
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.05+((figcount-1)*figdist) 0.45-((figcount-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.ple.fighandel(jk) = h;
end

