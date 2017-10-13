%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_bDT_einzel_RNAvsModel              %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Go through specified ligands and extract their prediction efficiency,
% F-measure and plot histograms and
% cumulative densitiy curves for the n training runs wrt random prediction

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

list_ligands = {'EGF','HRG','HGF','ALL'};%{'ALL'};%
BDT_plots = struct();
efficiency = NaN(2,length(list_ligands));

for i = 1:length(list_ligands)
    
    %Determine name of saved struct
    if(strcmp(list_ligands{i},'ALL'))
        lig_str = ['mut_wHGF_0_' num2str(perc)];
    else
        lig_str = [list_ligands{i} '_resp_0_' num2str(perc)];
    end
    
    %Calculate MPV = (PPV+NPV)/2 for both RNA and model based predictions
    efficiency(2,i) = plot_Conf(['RNA_' lig_str],save_folder,list_ligands{i});
    efficiency(1,i) = plot_Conf(['MODEL_' lig_str],save_folder,list_ligands{i});
    
    %Store MPV of all trained BDTs in BDT_plots struct
    BDT_plots.(['RNA_' list_ligands{i}]) = (bdt_figures.(['RNA_' lig_str])(1).testing(1).MPV(:,1,1) + bdt_figures.(['RNA_' lig_str])(1).testing(1).MPV(:,1,2))./2;
    BDT_plots.(['MODEL_' list_ligands{i}]) = (bdt_figures.(['MODEL_' lig_str])(1).testing(1).MPV(:,1,1) + bdt_figures.(['MODEL_' lig_str])(1).testing(1).MPV(:,1,2))./2;
    hold off
    for j=1:length(bdt_figures.(['MODEL_' lig_str]))
        %-0.5 to get it distributed around random prediction
        BDT_plots.(['MODEL_' list_ligands{i}])(j) = BDT_plots.(['MODEL_' list_ligands{i}])(j) - .5;
        BDT_plots.(['RNA_' list_ligands{i}])(j) = BDT_plots.(['RNA_' list_ligands{i}])(j) - .5;       
    end
    %Histogram plot of both MPV against random    
    figure
    histogram(BDT_plots.(['MODEL_' list_ligands{i}]),[-0.4:0.05:0.5])%,15)%
    hold on
    histogram(BDT_plots.(['RNA_' list_ligands{i}]),[-0.4:0.05:0.5])%,15)%
    set(gcf,'Color','w')
    xlabel('Difference BDT prediction to random [%]')
    ylabel('count')
    legend({'MODEL','RNA Expr'})
    hold off
    if(dosave)
        saveas(gcf,[save_folder 'BDT_hist_' list_ligands{i} '.fig'])
        saveas(gcf,[save_folder 'BDT_hist_' list_ligands{i}],'svg')        
    end
    %Cumulative density plot of both MPV against random
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
