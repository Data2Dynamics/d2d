%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   doBDT_calculations                      %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to calculate different BDT settings, model vs. RNAseq/qFACS
% with 1==model, 2==rna, and proportion of events taken for BDT testing

function doBDT_calculations(perc_testing,model_rna)
global bdt
global bdt_figures
nr_ligand = 4;

if(~exist('model_rna','var') || isempty(model_rna))
   model_rna = 1:2;
end

which_Structs = {'MODEL','RNA'};

% for i=1:length(which_Structs)
for i=model_rna
% del AB information
if(strcmp(which_Structs{i},'MODEL'))
    nr_entries = bdt.nr_model;
    do_model=1;
else
   nr_entries = bdt.nr_RNA; 
   do_model=0;
end
bdt.([which_Structs{i} '_Matrix'])(isnan(bdt.([which_Structs{i} '_Matrix'])(:,end)),:)=NaN;

% %delete IGF1 information
bdt.([which_Structs{i} '_Matrix'])((bdt.([which_Structs{i} '_Matrix'])(:,nr_entries+3)==1),:) = NaN;
bdt.([which_Structs{i} '_Matrix_bkp']) = bdt.([which_Structs{i} '_Matrix']);


bdt_figures.([which_Structs{i} '_EGF_resp_' strrep(num2str(perc_testing),'.','_')]) = BDT_bootstrap_cluster('EGF',do_model,500,perc_testing);
bdt_figures.([which_Structs{i} '_HRG_resp_' strrep(num2str(perc_testing),'.','_')]) = BDT_bootstrap_cluster('HRG',do_model,500,perc_testing);
bdt_figures.([which_Structs{i} '_HGF_resp_' strrep(num2str(perc_testing),'.','_')]) = BDT_bootstrap_cluster('HGF',do_model,500,perc_testing);

% bdt_figures.([which_Structs{i} '_IGF_resp']) = BDT_bootstrap_cluster('IGF1',do_model,500,0.4);
%continue;

% Cellular response to EGF, HRG, HGF without mutation information

% bdt.([which_Structs{i} '_Matrix'])(:,nr_entries-2:nr_entries) = NaN;
% bdt_figures.([which_Structs{i} '_nomut_wHGF']) = BDT_bootstrap_cluster('ALL',do_model,500,0.4);
% 
% %reset
% bdt.([which_Structs{i} '_Matrix']) = bdt.([which_Structs{i} '_Matrix_bkp']);

% continue;
% with mutation information, w/HGF

bdt_figures.([which_Structs{i} '_mut_wHGF_' strrep(num2str(perc_testing),'.','_')]) = BDT_bootstrap_cluster('ALL',do_model,500,perc_testing);

% %Delete everything except AKT/ERK
% if(do_model==1)
%     bdt.([which_Structs{i} '_Matrix'])(:,4:nr_entries-2) = NaN;
%     bdt_figures.([which_Structs{i} '_mut_woHGF_onlyDown']) = BDT_bootstrap_cluster('ALL',do_model,500,0.4);
% end
%reset
bdt.([which_Structs{i} '_Matrix']) = bdt.([which_Structs{i} '_Matrix_bkp']);

continue;
% % without PI3K cells
% bdt.([which_Structs{i} '_Matrix'])(bdt.([which_Structs{i} '_Matrix'])(:,nr_entries)==1,:) = NaN;
% bdt_figures.([which_Structs{i} '_noPI3K_wHGF']) = BDT_bootstrap_cluster('ALL',do_model,500,0.4);
% 
% %reset
% bdt.([which_Structs{i} '_Matrix']) = bdt.([which_Structs{i} '_Matrix_bkp']);
% 
% % without RAS cells
% bdt.([which_Structs{i} '_Matrix'])(bdt.([which_Structs{i} '_Matrix'])(:,nr_entries-1)==1,:) = NaN;
% bdt_figures.([which_Structs{i} '_noRAS_wHGF']) = BDT_bootstrap_cluster('ALL',do_model,500,0.4);
% 
% %reset
% bdt.([which_Structs{i} '_Matrix']) = bdt.([which_Structs{i} '_Matrix_bkp']);

end