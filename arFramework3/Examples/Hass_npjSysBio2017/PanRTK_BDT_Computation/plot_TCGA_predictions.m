%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_TCGA_predictions                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate all TCGA predictions, succeeded by barplots of predicted
% response of all indications to EGF, HRG or HGF

indications = {'BRCA','COADREAD','LUSC','LUAD','OV'};
ligands = {'HRG','EGF','HGF'};
%indications = {'PAAD'};
for i = 1:length(indications)
    for j = 1:length(ligands)
        name_tmp = [indications{i} '_' ligands{j}];
        bdt_name_tmp = ['TCGA_' indications{i}];
        if(strcmp(ligands{j},'HRG'))
            if(~run_regression)
                [TCGA_Matrices.(name_tmp), TCGA_Matrices.([name_tmp '_Prob'])] = predict(a,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+1)==1 & sum(bdt.(bdt_name_tmp)(:,bdt.nr_model+2:bdt.nr_model+4),2)==0,1:bdt.nr_model));            
            else
                TCGA_Matrices.([name_tmp '_Prob']) = predict(mdl_model,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+1)==1 & sum(bdt.(bdt_name_tmp)(:,bdt.nr_model+2:bdt.nr_model+4),2)==0,1:bdt.nr_model));
            end
        elseif(strcmp(ligands{j},'EGF'))            
            if(~run_regression)
                [TCGA_Matrices.(name_tmp), TCGA_Matrices.([name_tmp '_Prob'])] = predict(a,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+1)==0 & bdt.(bdt_name_tmp)(:,bdt.nr_model+4)==1,1:bdt.nr_model)); 
            else
                TCGA_Matrices.([name_tmp '_Prob']) = predict(mdl_model,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+1)==0 & bdt.(bdt_name_tmp)(:,bdt.nr_model+4)==1,1:bdt.nr_model));
            end
        elseif(strcmp(ligands{j},'HGF'))
            if(~run_regression)
                [TCGA_Matrices.(name_tmp), TCGA_Matrices.([name_tmp '_Prob'])] = predict(a,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+2)==1,1:bdt.nr_model));  
            else
                TCGA_Matrices.([name_tmp '_Prob']) = predict(mdl_model,bdt.(bdt_name_tmp)(bdt.(bdt_name_tmp)(:,bdt.nr_model+2)==1,1:bdt.nr_model));
            end
        end
        if(run_regression)
            TCGA_Matrices.(name_tmp) = double(TCGA_Matrices.([name_tmp '_Prob'])>opt_thresh);
        end
    end
end


save_folder='./TCGA_Predictions/';
% check if dir exists
if(~isempty(save_folder) && ~exist(save_folder, 'dir'))
    mkdir(save_folder)
end

%plot TCGA predictions separately for EGF/HRG/HGF
for j = 1:length(ligands)
    vec_tmp = [];
    for i = 1:length(indications)
        name_tmp = [indications{i} '_' ligands{j}];
        if(~run_regression)
            vec_tmp = [vec_tmp  sum(str2double(TCGA_Matrices.(name_tmp)))/length(TCGA_Matrices.(name_tmp))];
        else
            vec_tmp = [vec_tmp  sum(TCGA_Matrices.(name_tmp))/length(TCGA_Matrices.(name_tmp))];
        end
    end
    bar(1:length(indications),vec_tmp);
    title([ligands{j} ' stimulation'])
    h1 = get(gca,'Children');
    set(gca,'XLim',[0 length(indications)+1])
    set(gca,'XTick',1:length(indications),'XTickLabel',indications)
    set(gcf,'Color','w')
    ylabel('Cellular growth upon ligand stimulation [%]')
    saveas(gcf,[save_folder 'TCGA_' ligands{j} '_response.fig'])
    saveas(gcf,[save_folder 'TCGA_' ligands{j} '_response'],'svg')
end

