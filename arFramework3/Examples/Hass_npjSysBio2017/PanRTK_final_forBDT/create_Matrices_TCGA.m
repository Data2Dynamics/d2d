%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   create_Matrices_TCGA                    %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Load all TDGA indiciations, set offset to resemble expression level in
%  in vitro cell lines and calculate model response for every patient to be
%  fed into the BDT trained on the in vitro data
global TCGA_Matrices

uindications = unique(indications_select);

genes_list = {'EGFR','ERBB2','ERBB3','IGF1R','MET','EGF','NRG1','HGF','EGF','BTC','HBEGF','TGFA','AREG','EREG','EPGN'};

for i =1:length(uindications)
    var_tmp = ['TCGA_' uindications{i}];
    if(isfield(TCGA_Matrices,var_tmp))
        TCGA_Matrices.(var_tmp)=[];
    end
    for j=1:length(genes_list) 
        if(j<6)
            offset = mean(log10(sinh(data_select(strcmp(genes_select,genes_list{j}),:)))) - mean(log10(TCGA_Matrices.CCLEseq(j,:)));
            TCGA_Matrices.(var_tmp)(j,:) = 10.^(log10(sinh(data_select(strcmp(genes_select,genes_list{j}),ismember(indications_select,uindications{i}))))-offset);    
        else
            TCGA_Matrices.(var_tmp)(j,:) = sinh(data_select(strcmp(genes_select,genes_list{j}),ismember(indications_select,uindications{i})));
        end
    end
end

for i =1:length(uindications)
    ar.AB = [];
    var_tmp = ['TCGA_' uindications{i}];
    crawl_viability(TCGA_Matrices.(var_tmp)(1:5,:),1,1)
    set_RNA_Matrices(1,0,TCGA_Matrices.(var_tmp),var_tmp)
end