%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   plot_Conf                               %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting script for confusions matrix
% struct = the input struct in the global 'bdt_figures', i.e. the save name

% folder = output folder for the figures, if blank, BDT_plots is used

% which_lig = which ligand to be plotted, same as for BDT_bootstrap

function [hae, naive] = plot_Conf(struct, folder, which_lig) % [hae, naive] = 
global bdt_figures

if(~exist('struct','var') || isempty(struct))
    error('please provide the Classifier');
elseif(ischar(struct))
    struct_tmp = {struct};
else
    struct_tmp = struct;
end

if(~isfield(bdt_figures,struct))
    error('cannot find the field in bdt_figures struct');
end

if(~exist('folder','var') || isempty(folder))
    folder='BDT_plots';
end

% check if dir exists
if(folder~=0)
    if(~isempty(folder) && ~exist(folder, 'dir'))
        mkdir(folder)
    end
end
%find right struct for single MM/Ligand tests
which_drug = 'CONTROL';

if(~exist('which_lig','var') || isempty(which_lig) || (isinteger(which_lig) && which_lig==1))
    which_lig='ALL';
elseif((isinteger(which_lig) && which_lig==0))
    which_lig = 'CONTROL';
else
    which_lig = upper(which_lig );
end

if(~isempty(which_lig) && isempty(which_drug))
    which_drug='ALL';
elseif(~isempty(which_drug) && isempty(which_lig))
    which_lig='ALL';
end

if(~ischar(which_drug) || ~ischar(which_lig))
    error('please specify correct drug/ligand strings');
end

plot_id=[];

for i=1:length(struct_tmp)   
    if(~isempty(which_drug))
       for j=1:length(bdt_figures.(struct_tmp{i})(1).testing)
            if(strcmp(which_drug,bdt_figures.(struct_tmp{i})(1).testing(j).which_drugs) && strcmp(which_lig,bdt_figures.(struct_tmp{i})(1).testing(j).which_ligs))
               plot_id = j; 
            end
       end
    end
    if(isempty(plot_id))
        error('couldnt find the specified combination in the struct');
    end
    save_name1 = ['CONF_' struct_tmp{i}];
    save_name1 = ['./' folder '/' save_name1];
    if(~isempty(plot_id))
        save_name1= [save_name1 '_drug_' strrep(which_drug,'/',' and ') '_lig_' strrep(which_lig,'/',' and ') ];

    end
    titles={'Significant shrinkage','No sign. proliferation','Significant proliferation'};

    if(strfind(bdt_figures.(struct_tmp{i})(1).name,'RNA'))
        title_name = 'RNA, ';
    else
        title_name = 'model, ';
    end    
    if(~isempty(plot_id))
        title_name = [title_name 'drug ' strrep(which_drug,'_','\_') ', ligand ' which_lig];
    end
    %if(~isfield(bdt_figures.(struct_tmp{i})(1).testing(plot_id),'TFP') || isnan(bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(1,1,2)))
        
        %Store confusion matrix, normed by either number of events in
        %respective column or row
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed = NaN([3,3,length(bdt_figures.(struct_tmp{i}))] );
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed = NaN([3,3,length(bdt_figures.(struct_tmp{i}))] );
        
        %Store raw Conf Matrix
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_Matrix = NaN([3,3,length(bdt_figures.(struct_tmp{i}))] );
        
        %Store mean/std of normed Conf matrices
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean = NaN(3,3);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std = NaN(3,3);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean = NaN(3,3);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std = NaN(3,3);

        bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP = NaN(length(bdt_figures.(struct_tmp{i})),2,3);
        length_conf = size(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed,2);
    
        for j=1:length(bdt_figures.(struct_tmp{i}))
            %getting ugly, restruct cMat by hand
            C_mat_tmp = bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat;
            if(isempty(C_mat_tmp))
                continue;
            end
            Max_Cs = max(length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique),length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).bdt_class));
            if(Max_Cs>1)
                if(length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique)>length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).bdt_class))
                    missing_row = ismember(1:3,bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique+2);
                else
                    missing_row = ismember(1:3,bdt_figures.(struct_tmp{i})(j).testing(plot_id).bdt_class+2);
                end

                if(find(missing_row==0)==1)
                    C_mat_tmp = [zeros(1,3); zeros(2,1) bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat];
                elseif(find(missing_row==0)==3)
                    C_mat_tmp = [bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat zeros(2,1); zeros(1,3)];
                elseif(find(missing_row==0)==2)
                    C_mat_tmp = [bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat(1,1) 0 bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat(1,2); ...
                        zeros(1,3); bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat(2,1) 0 bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat(2,2)];
                end
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(:,:,j) = C_mat_tmp./(repmat(nansum(C_mat_tmp,1),length_conf,1));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(:,:,j) = C_mat_tmp./(repmat(nansum(C_mat_tmp,2),1,length_conf));
                
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_Matrix(:,:,j) = C_mat_tmp;
            elseif(Max_Cs==1)
                col_id = ismember([-1 0 1],bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique);
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(col_id,:,j) = 0;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(:,col_id,j) = 0;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(col_id,col_id,j) = 1;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(col_id,:,j) = 0;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(:,col_id,j) = 0;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(col_id,col_id,j) = 1;
            else
                continue;           
            end

        end
        for k=1:length_conf
            is_negative = 1:length_conf;
            is_negative(is_negative==k)=[];
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,k) = squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(k,k,:));
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,2,k) = squeeze(nansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(is_negative,k,:),1));
            for j=1:length_conf
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(j,k) = nanmean(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(j,k) = nanstd(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean(j,k) = nanmean(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std(j,k) = nanstd(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(j,k,:)));
            end

            mean_TP = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(k,k);
            std_TP = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(k,k);
            mean_FP = nansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(is_negative,k,:),1);
            std_FP = nansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(is_negative,k,:).^2,1);
            mean_FN = nansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean(k,is_negative,:),1);
            std_FN = nansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std(k,is_negative,:).^2,1);
            mean_recall = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean(k,k);
            std_recall = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std(k,k);

            ratio = mean_TP/(mean_FP);
            std_ratio = ratio * sqrt(std_TP^2/mean_TP^2 + std_FP/mean_FP^2);
            sprintf('Ratio of True/False for %s has mean %f +- %f',titles{k},ratio,std_ratio)
            sprintf('F-measure is: %f +- ',2*mean_TP*mean_recall/(mean_TP+mean_recall))

            sprintf('True predictions for it have mean %f +- %f',mean_TP,std_TP)
            sprintf('False predictions for it have mean %f +- %f',mean_FP,sqrt(std_FP))

            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).mean_TP = mean_TP;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_TP = std_TP;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).mean_FP = mean_FP;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_FP = std_FP;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).mean_FN = mean_FN;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_FN = std_FN;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).ratio = ratio;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_ratio = std_ratio;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_true = bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_1;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_false = sum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_0 + bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_neg1);
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_true = bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_0;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_false = sum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_1 + bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_neg1);
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_true = bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_neg1;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_false = sum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_0 + bdt_figures.(struct_tmp{i})(1).testing(plot_id).perc_1);

    %         [~, p_val] = ttest(squeeze(bdt_figures.(struct_tmp{i})(1).conf_col_normed(k,k,:)),squeeze(sum(bdt_figures.(struct_tmp{i})(1).conf_col_normed(is_negative,k,:),1)));
    %         sprintf('The t-test for %s has p-value of %f',titles{k},p_val)
        end
    %end
    hae =bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,3);
    naive = bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_true;
    
    sprintf('Averaged true predictions have mean %f +- %f',nanmean((bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,3) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,2))./2),nanstd((bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,3) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,2))./2))
    hae = nanmean((bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,3) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,2))./2);
    boxplot((bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,3) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP(:,1,2))./2,'labels',{'True predictions'})
    hold on
%     error_shrinkage = sqrt(bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_true*bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_false/bdt_figures.(struct_tmp{i})(1).testing(plot_id).events);
    error_neutral = sqrt(bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_true*bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_false/bdt_figures.(struct_tmp{i})(1).testing(plot_id).events); 
    error_growth = sqrt(bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_true*bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_false/bdt_figures.(struct_tmp{i})(1).testing(plot_id).events);    
    errorbar(1,.5,sqrt(error_growth^2+error_neutral^2),sqrt(error_growth^2+error_neutral^2),'Marker','.','Markersize',20,'Color','r')
    
%     %plot True/False assignment separately
%     boxplot(reshape(bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP,size(bdt_figures.(struct_tmp{i})(1).testing(plot_id).TFP,1),[],1),'labels',repmat({'True pred','False pred'},1,length_conf));    
%     hold on
%     Y_ax = get(gca,'YLim');
%     line([2.5 2.5],Y_ax)
%     line([4.5 4.5],Y_ax)
%     for j=1:length(titles)
%         text(1+(j-1)*2,Y_ax(1)+0.1,titles{j});
%     end
%     
%     errorbar(1,bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_true,error_shrinkage,error_shrinkage,'Marker','.','Markersize',20,'Color','k')
%     errorbar(2,bdt_figures.(struct_tmp{i})(1).testing(plot_id).neg_false,error_shrinkage,error_shrinkage,'Marker','.','Markersize',20,'Color','r')
%      
%     errorbar(3,bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_true,error_neutral,error_neutral,'Marker','.','Markersize',20,'Color','k')
%     errorbar(4,bdt_figures.(struct_tmp{i})(1).testing(plot_id).neutral_false,error_neutral,error_neutral,'Marker','.','Markersize',20,'Color','r')
%     
% %     errorbar(5,bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_true,error_growth,error_growth,'Marker','.','Markersize',20,'Color','k')
% %     errorbar(6,bdt_figures.(struct_tmp{i})(1).testing(plot_id).growth_false,error_growth,error_growth,'Marker','.','Markersize',20,'Color','r')
%     errorbar(5,.5,error_growth,error_growth,'Marker','.','Markersize',20,'Color','k')
%     errorbar(6,.5,error_growth,error_growth,'Marker','.','Markersize',20,'Color','r')

    h2 = get(gca,'Children');
%     legend([h2(1) h2(2)],'Binomial relative true data','Binomial relative false data')
%     ylabel('Relative true/false predictions','Interpreter','Tex');
    legend([h2(1)],'Binomial relative true data')
    ylabel('Relative true predictions','Interpreter','Tex');
    
    set(gcf,'Color','w');
    title(sprintf('Statistics of training for %s',title_name));
    if(folder~=0)
        saveas(gcf,save_name1,'fig')
        saveas(gcf,save_name1,'pdf')
    end
end

