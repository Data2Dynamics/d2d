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
% MPV = mean predictive value, the measure we look at (PPV+NPV)/2
% The matrix in this script has MPV = [NPV_BDT1 PPV_BDT1]
%                                     [NPV_BDT2 PPV_BDT2]
%                                     [...      ...     ]
%
% ConfusionMatrix as given from BDT_getROC (0 = no growth, 1 = growth)
% ROWS=TRUE VALUES / COLS= PRED VALUES
%    0    1
% 0 TN   FP
% 1 FN   TP
%
% For each predicted value, the following measures are calculated as mean of all BDTs and
% stored in the bdt_figures struct Nr. 1:
%  mean_PV = positive/negative predictive value
%  mean_FDR = false discovery rate
%  mean_recall = recall


function MPV = plot_Conf(struct, folder, which_lig) % [hae, naive] = 
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
%find right struct for single Ligand tests
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
    titles={'No sign. proliferation','Significant proliferation'};

    if(strfind(bdt_figures.(struct_tmp{i})(1).name,'RNA'))
        title_name = 'RNA, ';
    else
        title_name = 'model, ';
    end    
    if(~isempty(plot_id))
        title_name = [title_name 'drug ' strrep(which_drug,'_','\_') ', ligand ' which_lig];
    end
    %if(~isfield(bdt_figures.(struct_tmp{i})(1).testing(plot_id),'MPV') || isnan(bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(1,1,2)))
        
        %Store confusion matrix, normed by either number of events in
        %respective column or row
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed = NaN([2,2,length(bdt_figures.(struct_tmp{i}))] );
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed = NaN([2,2,length(bdt_figures.(struct_tmp{i}))] );
        
        %Store raw Conf Matrix
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_Matrix = NaN([2,2,length(bdt_figures.(struct_tmp{i}))] );
        
        %Store mean/std over all trained BDTs of normed Conf matrices
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean = NaN(2,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std = NaN(2,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean = NaN(2,2);
        bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std = NaN(2,2);

        bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV = NaN(length(bdt_figures.(struct_tmp{i})),1,2);
        length_conf = size(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed,2);
        for j=1:length(bdt_figures.(struct_tmp{i}))
        %loop overall trained BDTs, and store Matrices normed by sum over
        %respectrive row (col_normed) or column (row_normed)
            C_mat_tmp = bdt_figures.(struct_tmp{i})(j).testing(plot_id).C_mat;
            if(isempty(C_mat_tmp))
                continue;
            end
            Max_Cs = max(length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).Y_unique),length(bdt_figures.(struct_tmp{i})(j).testing(plot_id).bdt_class));
            if(Max_Cs>1)                
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(:,:,j) = C_mat_tmp./(repmat(arnansum(C_mat_tmp,1),length_conf,1));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(:,:,j) = C_mat_tmp./(repmat(arnansum(C_mat_tmp,2),1,length_conf));
                
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_Matrix(:,:,j) = C_mat_tmp;
            elseif(Max_Cs==1)
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(:) = NaN;
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(:) = NaN;
            else
                continue;           
            end

        end
        for k=1:length_conf
            %This loop is needed if more than 0/1 decisions have to be
            %made. Here, means over all BDTs are calculated for the normed
            %confusion matrices stored above
            is_negative = 1:length_conf;
            is_negative(is_negative==k)=[];
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,k) = squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(k,k,:));
            for j=1:length_conf
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(j,k) = nanmean(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(j,k) = nanstd(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean(j,k) = nanmean(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(j,k,:)));
                bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std(j,k) = nanstd(squeeze(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_normed(j,k,:)));
            end
            %Some statistical properties are calculated (see description in
            %the header)
            mean_PV = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(k,k);
            std_PV = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(k,k);
            mean_FDR = arnansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_mean(is_negative,k,:),1);
            std_FDR = arnansum(bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_col_std(is_negative,k,:).^2,1);
            mean_recall = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_mean(k,k);
            std_recall = bdt_figures.(struct_tmp{i})(1).testing(plot_id).conf_row_std(k,k);

            ratio = mean_PV/(mean_FDR);
            std_ratio = ratio * sqrt(std_PV^2/mean_PV^2 + std_FDR/mean_FDR^2);
            sprintf('Ratio of True/False for %s has mean %f +- %f',titles{k},ratio,std_ratio)
            sprintf('F-measure is: %f +- ',2*mean_PV*mean_recall/(mean_PV+mean_recall))

            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).mean_TP = mean_PV;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_TP = std_PV;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).mean_FP = mean_FDR;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_FP = std_FDR;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).ratio = ratio;
            bdt_figures.(struct_tmp{i})(1).testing(plot_id).Conf(k).std_ratio = std_ratio;
        end
    %Function returns MPV = (PPV+NPV)/2 and prints it
    MPV = nanmean((bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,1) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,2))./2);    
    sprintf('Averaged true predictions have mean %f +- %f',nanmean((bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,1) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,2))./2),nanstd((bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,1) + bdt_figures.(struct_tmp{i})(1).testing(plot_id).MPV(:,1,2))./2))

end

