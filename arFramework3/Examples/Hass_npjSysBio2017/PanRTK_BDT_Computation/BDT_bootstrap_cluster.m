%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   BDT_bootstrap_cluster                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to run BDTs on viability data for single or multiple ligands/ABs
% with parallel computing
% whichLig = 'ALL' to test on all ligands, 'ITER' to loop through all
% ligands and test the single ligand response, or e.g. 'HRG', or 'HRG/IGF1'
% to get response from one specific ligand.

% for both, specify 'control' to get the BDT to the control condition
% e.g. whichLig='all', whichMM='control' to get ligand effect of
% combinations of all ligands

% ligand names = 'HRG','HGF','IGF1','EGF'

% model=0 specifies RNA, model=1 if model output should be used
% nboot is number of trees trained
% perc_testing is the percentage of cell lines used for testing

% do_rnd = 1 trains BDT on randomized output, should give 50:50

function classifier = BDT_bootstrap_cluster(whichLig, model, nboot, perc_testing, do_rnd)
global bdt
if(~exist('nboot','var') || isempty(nboot))
    nboot=500;
end
if(~exist('perc_testing','var') || isempty(perc_testing))
    perc_testing = 0.2;
end
if(~exist('do_rnd','var') || isempty(do_rnd))
    do_rnd=0;
end

if(~exist('whichLig','var') || isempty(whichLig) || (isnumeric(whichLig) && whichLig==0))
    whichLig='ALL';
elseif((isnumeric(whichLig) && whichLig==1))
    whichLig = 'ITER';
elseif(ischar(whichLig))
    whichLig = upper(whichLig);
else
    error('please specifiy Ligand either by 0/1, empty or string');
end

which{1} = 'RNA_Matrix';
if(model)
    which{1}='MODEL_Matrix';
end

%This sets equal cell line samples in each run
%Comment out the respective lines in the beginning of parfor loop
spmd
    rng(16,'combRecursive');
end

%How many minimum leafs per training (3), number of trees consecutively
%trained (300)
leaf=5;
nTrees = 300;

%Predictors X, Response Y
X = bdt.(which{1})(:,1:end-1);
Y = bdt.(which{1})(:,end);

%How many cell lines are there?
nCellLines = length(bdt.cell_names);
%get repetitions of one cell line in RNA sample
nr_reps=bdt.nr_reps;
xVal = [];
yVal = [];
auc = [];
startTime = clock;
arShowProgressParFor(nboot);
classifier=struct();
bdt1 = bdt;
if(~model)
    nr_RNA = bdt.nr_RNA;
else
    nr_RNA = bdt.nr_model;
end
    
parfor i=1:nboot
    bdt_tmp = bdt1;
    
    %Different random setups. Comment out to set split of
    %training/testing equal in every run of this script
    stream = RandStream.getGlobalStream();
    stream.Substream = i;

    classifier(i).name = which{1};
    
    %Determine which entries are used for training/testing
    %cells_out gives array of the cell lines used for testing
    cells_out_IDs = randi([1 nCellLines],[1 floor(nCellLines*perc_testing)]);

    add_matrix = repmat([1:nr_reps]',1,length(cells_out_IDs));         
    
    cells_out = repmat((cells_out_IDs-1)*size(add_matrix,1),size(add_matrix,1),1) + add_matrix;
    cells_out = reshape(cells_out,[],1);

    %ugly procedure to get training/testing/random sets, delete all nan's
    X_training = X;
    X_training(cells_out,:) = [];
    Y_training = Y;
    Y_training(cells_out) = [];
    
    if(~model)
        nr_RNA = bdt_tmp.nr_RNA;
    else
        nr_RNA = bdt_tmp.nr_model;
    end
    nr_ligs = length(bdt_tmp.lig_names);

    [X_training, Y_training] = del_nans(X_training, Y_training);
    if(model)
        X_training = X_training(:,1:nr_RNA);
    end
    
     X_test = X(cells_out,:);
     Y_test = Y(cells_out,:);
    [X_test, Y_test] = del_nans(X_test, Y_test);
    b_tmp = TreeBagger(nTrees, X_training, Y_training, 'SampleWithReplacement','on',...
           'Method','classification','MinLeafSize',leaf,'FBoot',0.8);%,'UseParallel',true); 
    
%     classifier(i).tree = b_tmp;
    classifier(i).tree_lig = whichLig;
    classifier(i).tree_drug = 'CONTROL';
    
    %do it for the randomized sample
    if(do_rnd)
        %Get random data        
        shuffle_Y = randperm(length(Y_training));
        Y_rnd_training = Y_training(shuffle_Y);        
        X_rnd_training = X_training;

        b_tmp_rnd = TreeBagger(nTrees, X_rnd_training, Y_rnd_training, 'SampleWithReplacement','on',...
           'Method','classification','MinLeafSize',leaf,'FBoot',0.8);%,'UseParallel',true);        
        classifier(i).tree_rnd = b_tmp_rnd;
    end
    
    %puzzle together Ligand combinations
    
    if(~strcmp(whichLig,'ALL') && ~strcmp(whichLig,'ITER'))
        which_ligs = get_input_string(bdt_tmp,whichLig);
    elseif(strcmp(whichLig,'ITER'))
        which_ligs = unique(X_test(:,nr_RNA+1:nr_RNA+nr_ligs),'rows');
    else
        which_ligs = 0;
    end
    isave=1;
    for k=1:size(which_ligs,1)
        isave = isave+1;
        lig_tmp = which_ligs(k,:);

        %find correct names of ligands
        if(length(lig_tmp)==1)
            lig_names='ALL';
        else
            lig_names = bdt_tmp.lig_names(find(lig_tmp==1));
            if(isempty(lig_names))
                lig_names = 'CONTROL'; 
            else
                lig_names = strjoin(lig_names,'/');
            end
        end

        [X_test_tmp, Y_test_tmp] = select_test(model,bdt_tmp, X_test, Y_test, lig_tmp);
        [xVal, yVal, auc, confMat, bdt_classes, Class_prob] = BDT_getROC(b_tmp,X_test_tmp,Y_test_tmp);

        if(i==1)
            [X_trivial, Y_trivial] = del_nans(X,Y);                  
            [~, Y_trivial] = select_test(model,bdt_tmp,X_trivial,Y_trivial,lig_tmp);
            classifier(i).testing(isave).perc_1 = sum(Y_trivial==1)/length(Y_trivial);
            classifier(i).testing(isave).perc_0 = sum(Y_trivial==0)/length(Y_trivial);
            classifier(i).testing(isave).perc_neg1 = sum(Y_trivial==-1)/length(Y_trivial);
            classifier(i).testing(isave).events = length(Y_trivial);
        end
        classifier(i).testing(isave).xVal = xVal;
        classifier(i).testing(isave).yVal = yVal;
        classifier(i).testing(isave).auc = auc;
        classifier(i).testing(isave).C_mat = confMat;            
        classifier(i).testing(isave).which_ligs = lig_names;
        classifier(i).testing(isave).which_drugs = 'CONTROL';
        classifier(i).testing(isave).bdt_class = bdt_classes;
        classifier(i).testing(isave).Y_unique = unique(Y_test_tmp);
        classifier(i).naive_growth = sum(Y_test_tmp==1)/length(Y_test_tmp);

        if(do_rnd)
            Y_rnd_test = Y_test_tmp(randperm(length(Y_test_tmp)));
            X_rnd_test=X_test_tmp;        
           [xVal, yVal, auc, confMat, bdt_rnd_classes, ~] = BDT_getROC(b_tmp_rnd,X_rnd_test,Y_rnd_test);
           %bdt_figures.(saveName)(i).Class_rnd=b_tmp;
           classifier(i).testing(isave).C_rnd_mat = confMat;   
           classifier(i).testing(isave).xVal_rnd = xVal;
           classifier(i).testing(isave).yVal_rnd = yVal;
           classifier(i).testing(isave).auc_rnd = auc;
           classifier(i).testing(isave).bdt_rnd_class = bdt_rnd_classes;
           classifier(i).testing(isave).Y_rnd_unique = unique(Y_rnd_test);
        end
    end
       
  
   
  
    return_message = sprintf('im in loop %i', i);
    arShowProgressParFor(i, nboot, startTime, return_message)
end
arShowProgressParFor(0);


function [out_X, out_Y] = del_nans(in_X, in_Y)
    out_X = in_X;
    out_Y = in_Y;
    out_X(isnan(out_Y),:)=[];
    out_Y(isnan(out_Y))=[];
    out_Y(sum(isnan(out_X),2)==size(out_X,2))=[];
    out_X(sum(isnan(out_X),2)==size(out_X,2),:)=[];


function [out_X, out_Y] = select_test(model,bdt_tmp, in_X,in_Y, which_ligs)

if(~model)
    nr_RNA = bdt_tmp.nr_RNA;
else
    nr_RNA = bdt_tmp.nr_model;
end
nr_ligs = length(bdt_tmp.lig_names);

if(length(which_ligs) == 1)
    out_X = in_X;
    out_Y = in_Y;
else
    out_X = in_X(ismember(in_X(:,nr_RNA+1:nr_RNA+nr_ligs),which_ligs,'rows'),:);
    out_Y = in_Y(ismember(in_X(:,nr_RNA+1:nr_RNA+nr_ligs),which_ligs,'rows'));
end
if(length(out_Y)<10)
   sprintf('There are just %d events or extremely low event count for your specified condition. \n',length(out_Y))
   sprintf('Training/testing will continue since you can re-use the trained trees. \n')      
end
if(model)
    out_X = out_X(:,1:nr_RNA);
end

function which_lig = get_input_string(bdt_tmp,input)
%get input strings

    struct_name = 'lig_names';
    
    if(ischar(input))
       double = strsplit(input,'/');
       if(length(double)==1)
          which_lig = strcmp(bdt_tmp.(struct_name),double{1})+0;
       elseif(length(double)==2)
           if(strcmp(double{1},double{2}))
               error('Haha');
           end
          which_lig = zeros(1,length(bdt_tmp.(struct_name)));
          for nd=1:length(double)
             which_lig = which_lig + strcmp(bdt_tmp.(struct_name),double{nd})+0;
          end
       else
           error('only single or double lig combos!');
       end
       if(sum(which_lig)==0 && ~strcmp(input,'CONTROL'))
          error('Please specify existing ligands or ALL, ITER, CONTROL'); 
       end
    end