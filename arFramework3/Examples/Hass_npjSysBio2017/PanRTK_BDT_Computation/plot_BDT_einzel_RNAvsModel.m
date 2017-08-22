%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_bDT_einzel_RNAvsModel              %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Go through specified ligands and extract their prediction efficiency,
% F-measure and plot boxplots of TP/FP ratio as well as histograms and
% cumulative densitiy curves for the n training runs wrt naive growth

function efficiency = plot_BDT_einzel_RNAvsModel(dosave,perc)
global bdt_figures
global BDT_plots
save_folder =0;
if(~exist('perc','var') || isempty(perc))
   perc = 3;
end
if(dosave)
    save_folder='./BDT_plots/';
    % check if dir exists
    if(~isempty(save_folder) && ~exist(save_folder, 'dir'))
        mkdir(save_folder)
    end
end

list_ligands = {'EGF','HRG','HGF','ALL'}; %'IGF1','HRG/IGF1'
BDT_plots = struct();
efficiency = NaN(2,length(list_ligands));

for i = 1:length(list_ligands)
    
    if(strcmp(list_ligands{i},'ALL'))
        lig_str = ['mut_wHGF_0_' num2str(perc)];
    else
        lig_str = [list_ligands{i} '_resp_0_' num2str(perc)];
    end
    figure
    efficiency(2,i) = plot_Conf(['RNA_' lig_str],save_folder,list_ligands{i});
    figure
    efficiency(1,i) = plot_Conf(['MODEL_' lig_str],save_folder,list_ligands{i});
%     BDT_plots.(['RNA_' list_ligands{i}]) = bdt_figures.(['RNA_' lig_str])(1).testing(2).TFP(:,1,3);
%     BDT_plots.(['MODEL_' list_ligands{i}]) = bdt_figures.(['MODEL_' lig_str])(1).testing(2).TFP(:,1,3);
    BDT_plots.(['RNA_' list_ligands{i}]) = (bdt_figures.(['RNA_' lig_str])(1).testing(2).TFP(:,1,3) + bdt_figures.(['RNA_' lig_str])(1).testing(2).TFP(:,1,2))./2;
    BDT_plots.(['MODEL_' list_ligands{i}]) = (bdt_figures.(['MODEL_' lig_str])(1).testing(2).TFP(:,1,3) + bdt_figures.(['MODEL_' lig_str])(1).testing(2).TFP(:,1,2))./2;
    hold off
    for j=1:length(bdt_figures.(['MODEL_' lig_str]))
%         BDT_plots.(['MODEL_' list_ligands{i}])(j) = BDT_plots.(['MODEL_' list_ligands{i}])(j) - bdt_figures.(['MODEL_' lig_str])(j).naive_growth;
%         BDT_plots.(['RNA_' list_ligands{i}])(j) = BDT_plots.(['RNA_' list_ligands{i}])(j) - bdt_figures.(['RNA_' lig_str])(j).naive_growth;       
        BDT_plots.(['MODEL_' list_ligands{i}])(j) = BDT_plots.(['MODEL_' list_ligands{i}])(j) - .5;
        BDT_plots.(['RNA_' list_ligands{i}])(j) = BDT_plots.(['RNA_' list_ligands{i}])(j) - .5;       
    end
        
    figure
    histogram(BDT_plots.(['MODEL_' list_ligands{i}]),15)
    hold on
    histogram(BDT_plots.(['RNA_' list_ligands{i}]),15)
    set(gcf,'Color','w')
    % set(gca,'YLim',[0 100])
    xlabel('Difference BDT prediction to random [%]')
    ylabel('count')
    legend({'MODEL','RNA Expr'})
    hold off
    if(dosave)
        saveas(gcf,[save_folder 'BDT_hist_' list_ligands{i} '.fig'])
        saveas(gcf,[save_folder 'BDT_hist_' list_ligands{i}],'svg')        
    end
    
    figure
    cdfplot(BDT_plots.(['MODEL_' list_ligands{i}]))
    hold on
    cdfplot(BDT_plots.(['RNA_' list_ligands{i}]))
    set(gcf,'Color','w')
    legend({'MODEL','RNA Expr'})
    xlabel('Difference BDT prediction to random [%]')
    ylabel('Cumulative distribution')
    title('')
    hold off
    if(dosave)
        saveas(gcf,[save_folder 'BDT_cdf_' list_ligands{i} '.fig'])
        saveas(gcf,[save_folder 'BDT_cdf_' list_ligands{i}],'svg')        
    end
end
