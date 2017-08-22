%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   set_RNA_Matrices                        %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  script to create matrices with RNAseq and binary ligand columns or model
%  features mutations are added as well and viability outcome is added if
%  the matrices are created for in-vitro cell lines


function set_RNA_Matrices(model,doY,RNA,Matrix_name)

global Matrices
global bdt
global bdt_figures
global ar

if(~exist('doY','var') || isempty(doY))
    doY = 1;
end

if(~exist('RNA','var') || isempty(RNA))
    %For all cell lines, repeat RNAseq N times, with N coming from the size of
    %viability matrix (Control vs. 6 ligands)
    %nr_cell2s has rows: EGFR, ErbB2, ErbB3, IGF1R, ERBB4, MET, HGF, EGF, HRG, IGF
%     x_RNA_tmp = sqrt(Matrices.RNAseq([1:4 6],:)');
    x_RNA_tmp = (Matrices.qFACS(1:5,:)/10000)';
%     x_RNA_tmp = [(Matrices.qFACS(1:5,:)/10000)' Matrices.RNAseq([8 9 7 10],:)'];
else
    x_RNA_tmp = RNA';
end
if(~exist('Matrix_name','var') || isempty(Matrix_name))
    if(~model)
        Matrix_name = 'RNA_Matrix';
    else
        Matrix_name = 'MODEL_Matrix';
    end
end

bdt.cell_names = Matrices.RNAseq_cellnames;

bdt.lig_names = {'HRG','HGF','IGF1','EGF'};
%mutation columns: KRAS, PI3K
mut_cols = 2;

nligs = size(Matrices.Dflat_vsLig,2);
ncells = size(x_RNA_tmp,1);
%ncells = 1000;
% x_RNA_LigvsCtrl = bdt.MODEL_Matrix(:,1:end-1);
nRNA = size(x_RNA_tmp,2);

bdt.nr_RNA = nRNA+mut_cols;

list_KRAS = {'GP2D','H747','HCT116','HCT15','LOVO','LS123','LS180','RCM-1','SW620','T84','H441','A549','H358','MDA-MB-231','AGS','H23','H460'};
list_pi3k = {'CCK81','GP2D','H508','HCT116','HCT15','HT115','LS180','RKO','SW48','T84','BT474','HCC1954','T47D','AGS','H460','H596'};

%create Matrix for Ligand effect only
x_RNA_LigvsCtrl = [];   

%Set RNAseq in matrices
if(~model)
    for i=1:ncells
       x_RNA_LigvsCtrl = [x_RNA_LigvsCtrl; repmat(x_RNA_tmp(i,:),nligs-1,1)]; %nligs-1 since 'control' is not taken in matrix
    end
else
    model_input = fetch_outputs;
    
    %do the same for model output
%Set Model output in matrices
    for i=1:length(model_input)
        x_RNA_LigvsCtrl =  [x_RNA_LigvsCtrl reshape(ar.AB.(model_input{i})(1,2:end,:),[],1,1)];
    end
    nligs = size(ar.AB.pS6,2);
    nRNA = length(model_input);

    bdt.nr_model = nRNA+mut_cols;
end
length_LigvsCtrl = ncells*(nligs-1);
%Expand matrices by 4 columns for HRG, HGF, IGF1, EGF

x_RNA_LigvsCtrl = [x_RNA_LigvsCtrl zeros(length_LigvsCtrl,mut_cols) zeros(length_LigvsCtrl,4)];

%Viability structure
LigvsCtrl_rows_perCell = nligs-1;

if(~exist('RNA','var') || isempty(RNA))
    %fill mutations status
    has_PI3K = find(ismember(bdt.cell_names,list_pi3k));
    has_RAS = find(ismember(bdt.cell_names,list_KRAS));    
    
    x_RNA_LigvsCtrl(reshape(repmat((has_RAS-1)*LigvsCtrl_rows_perCell,1,LigvsCtrl_rows_perCell)+repmat(1:LigvsCtrl_rows_perCell,size(has_RAS,1),1),[],1),nRNA+1)=1;
    x_RNA_LigvsCtrl(reshape(repmat((has_PI3K-1)*LigvsCtrl_rows_perCell,1,LigvsCtrl_rows_perCell)+repmat(1:LigvsCtrl_rows_perCell,size(has_PI3K,1),1),[],1),nRNA+2)=1;
end

%Mapping Dflat structure to ligand columns in Prediction-matrix
ligs_cols = [nRNA+mut_cols+1 NaN; ...
    nRNA+mut_cols+2 NaN; ...
    nRNA+mut_cols+3 NaN; ...
    nRNA+mut_cols+1 nRNA+mut_cols+3; ... %HRG/IGF1
    nRNA+mut_cols+4 NaN; ...
    nRNA+mut_cols+1 nRNA+mut_cols+4];  %HRG/EGF

%set Ligand columns
for i = 1:nligs-1 %nligs-1 since first row is skipped (no ligand applied)
    x_RNA_LigvsCtrl(i:LigvsCtrl_rows_perCell:length_LigvsCtrl,ligs_cols(i,~isnan(ligs_cols(i,:)))) = 1;    
end

if(doY)
    %get column of outputs, set categories for 20% growth/shrinkage and append
    %save output matrices in Matrices struct for sanity checks
   
    Matrices.Dsig_LigvsCtrl_Matrix = NaN(1,nligs-1,ncells);
    Matrices.Dflat_LigvsCtrl_Matrix = NaN(1,nligs-1,ncells);

    for j=1:length(Matrices.RNAseq_cellnames)
           for i=1:length(Matrices.Exp_names)
                if(strcmp(Matrices.RNAseq_cellnames{j},Matrices.Exp_names{i}))                   
                    Matrices.Dsig_LigvsCtrl_Matrix(:,:,j) = Matrices.Dsig_vsCtrl(1,2:end,i); %_10
                    Matrices.Dflat_LigvsCtrl_Matrix(:,:,j) = Matrices.Dflat_vsCtrl(1,2:end,i)-1; %_10
                end
           end
    end

    Dflats = {'Dflat_LigvsCtrl_Matrix'};
    Dsigs = {'Dsig_LigvsCtrl_Matrix'};

    if(~model)
        name_ext = 'x_RNA_';
    else
        name_ext = 'x_model_';
    end

    for i=1:length(Dflats)
        D_name = strsplit(Dflats{i},'_');
        D_name = D_name{2};
        y_tmp = reshape(Matrices.(Dflats{i}),[],1,1);
        sig_tmp = reshape(Matrices.(Dsigs{i}),[],1,1);

        y_RNA_tmp = y_tmp;
        y_RNA_tmp(~isnan(y_tmp)) = 0;
        y_RNA_tmp(sig_tmp==1) = 1;
        y_RNA_tmp(y_tmp<0 & sig_tmp==1) = -1;

        x_name = ['x_RNA_' D_name];
        new_name = [name_ext D_name];

        matrix_tmp = [eval(x_name) y_RNA_tmp];
        Matrices.(new_name) = matrix_tmp;

    end
end

bdt.nr_reps = LigvsCtrl_rows_perCell;
%Set combination matrix with all data
if(doY)
    bdt.(Matrix_name) = Matrices.([name_ext 'LigvsCtrl']);
else
    bdt.(Matrix_name) = x_RNA_LigvsCtrl;
end

if(~doY)
    bdt.(Matrix_name)(isnan(bdt.(Matrix_name)(:,1)),:) = [];
end

