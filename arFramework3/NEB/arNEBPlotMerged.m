function arNEBPlotMerged
%Plots result of the merging approach

global ar

subplot(2,1,1)

chi2_lhs_plot = (ar.chi2s_sorted-(ar.chi2s_sorted(1) -1));
chi2_start = (ar.chi2s_start_sorted-(ar.chi2s_sorted(1)) -1);

semilogy(chi2_start,'--x', 'Color', [0.6, 0.6, 0.6] )
hold on
semilogy(chi2_lhs_plot,'-o', 'Color', [1, 0, 0 ], 'Markersize',8 , 'MarkerFaceColor', [1, 0, 0 ])
% semilogy(chi2_lhs_plot,'-o', 'Color', [1, 0.3, 0 ], 'Markersize',8 , 'MarkerFaceColor', [1, 0.3, 0 ])
hold on


con_min = chi2_lhs_plot(ar.merger.merged_minima_neb);

springs = unique(ar.merger.merged_minima_spring_neb);
springs = springs(~isnan(springs));
springs = springs(~(springs==0));

col_spr = winter(length(springs));

for i = 1:length(chi2_lhs_plot)
   
    if isnan(ar.merger.merged_minima_spring_neb(i))
        
        semilogy(i, chi2_lhs_plot(ar.merger.merged_minima_neb(i)),'xk', 'Markersize',12 )

    elseif ar.merger.merged_minima_spring_neb(i) == 0
            
        semilogy(i, chi2_lhs_plot(ar.merger.merged_minima_neb(i)),'ok', 'Markersize',8 )

    else
        [~,xx] = find(springs == ar.merger.merged_minima_spring_neb(i) );
        semilogy(i, chi2_lhs_plot(ar.merger.merged_minima_neb(i)),'ok', 'Color', col_spr(xx,:), 'Markersize',8 , 'MarkerFaceColor', col_spr(xx,:))      

    end
    
    
end



for i = 1:length(con_min)
    line( [i i], [con_min(i) chi2_lhs_plot(i)], 'Color', [0, 0, 0.8], 'LineWidth',2 )
    line( [ar.merger.merged_minima_neb(i) i] , [con_min(i) con_min(i)], 'Color', [0.5, 0.5, 0.5], 'LineWidth',1 )
end

%uni_con_min = unique(con_min);

ylims = ylim;

ylim([0.4 ylims(2)])

plot([1, 50], [1,1], '--k')


ylabel('likelihood')
xlabel('fit index sorted by (suboptimal) multistart optimization result')
title('before merging')

legend({'initial guess','after fitting (regular LHS result)','after merging ','connection to minimum via optimal path ',})



subplot(2,1,2)

[x,xx] = sort(chi2_lhs_plot(ar.merger.merged_minima_neb));
semilogy(chi2_start(xx),'--x', 'Color', [0.6, 0.6, 0.6] )
hold on

semilogy(chi2_lhs_plot(xx),'o', 'Color', [1, 0, 0 ], 'Markersize',5, 'MarkerFaceColor', [1, 0, 0 ])

ylims = ylim;

semilogy(x, 'o-', 'Color', [0, 0, 0.8], 'Markersize',6 , 'MarkerFaceColor', [0, 0, 0.8] )


[y,yy] = unique(x);

for i = 1:length(y)
     %line([yy(i) 50],[y(i) y(i) ], 'Color', [0.4, 0.4, 0.4])
     %plot( [yy(i) 50]', [y(i) y(i) ]', 'Color', [0.4, 0.4, 0.4],'--' )
     plot([yy(i), 50], [y(i),y(i)], '-k')

end

ylim([0.4 ylims(2)])

ylabel('likelihood')
xlabel('fit index sorted after merging of connected minima')
title('after merging')

legend({'initial guess','after fitting (regular LHS result)','after merging ',' local minima',})

end





