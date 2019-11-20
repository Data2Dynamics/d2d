function arNEBPlotpair
% Plots pairwise paths in parameter space for the current NEB for all
% spring constants

global ar

figure

ncol = length(ar.merger.neb.qFitlable);
nrow = length(ar.merger.neb.qFitlable);

legcell = {};

for i_spr = 1:length(ar.merger.neb.spr)
    
    
    for i_row = 1:nrow
        for i_col = 1:ncol
            
            if i_row == 1
                subplot(nrow,ncol,(i_row-1)*ncol+i_col)
                plot(ar.merger.neb.spr(i_spr).ps_result(:,i_col),'-')
                hold on
                title(ar.merger.neb.qFitlable(i_col))
            else
                if i_col < i_row
                	subplot(nrow,ncol,(i_row-1)*ncol+i_col)
                    plot(ar.merger.neb.spr(i_spr).ps_result(:,i_col),ar.merger.neb.spr(i_spr).ps_result(:,i_row),'-')
                    hold on
                    
                    plot(ar.merger.neb.spr(i_spr).ps_result(1,i_col),ar.merger.neb.spr(i_spr).ps_result(1,i_row),'xk')
                    plot(ar.merger.neb.spr(i_spr).ps_result(end,i_col),ar.merger.neb.spr(i_spr).ps_result(end,i_row),'xk')
                    if i_col == 1
                        ylabel(ar.merger.neb.qFitlable(i_row))
                    end
                    
                    if i_row == ncol
                        xlabel(ar.merger.neb.qFitlable(i_col))
                    end 
                end     
            end
        end
    end
    
    legcell{i_spr} =  ['k = ' num2str(ar.merger.neb.spr(i_spr).springconst)];
    subplot(nrow,ncol,3*ncol)
    plot( 1, ar.merger.neb.spr(i_spr).springconst,'-o' )
    hold on
    legend(legcell,'Location','southoutside')
end

    
end
