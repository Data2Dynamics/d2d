%small script to train BDT based on all in-vitro events and plot 2D slices
%of features space
addpath('../PanRTK_final_forBDT')
list_outputs = fetch_outputs;
D_list = {'pEGFR_hom','pEGFR_12','pS6','pAKT','pERK','pMet_ErbB3'};

Mat_MODEL_mut = bdt.MODEL_Matrix;
Mat_MODEL_mut(isnan(Mat_MODEL_mut(:,end)),:)=[];

id_true = find(Mat_MODEL_mut(:,end)==1);
id_false = find(Mat_MODEL_mut(:,end)==0);

scatter(log(Mat_MODEL_mut(id_false,find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_false,find(ismember(list_outputs,D_list{2})))),(1.05.^(Mat_MODEL_mut(id_false,1))+0)*12,zeros(length(id_false),find(ismember(list_outputs,D_list{3}))),'filled')
hold on
scatter(log(Mat_MODEL_mut(id_true,find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_true,find(ismember(list_outputs,D_list{2})))),(1.05.^(Mat_MODEL_mut(id_true,1))+0)*12,ones(length(id_true),find(ismember(list_outputs,D_list{3}))),'filled')
set(gcf,'Color','w')
colormap(jet)
xlabel('EGFR homodimerization')
ylabel('EGFR-ErbB2 heterodimerization')

id_highErbB2 = find(log(Mat_MODEL_mut(:,find(ismember(list_outputs,D_list{2}))))>3);
id_lowErbB2 = find(log(Mat_MODEL_mut(:,find(ismember(list_outputs,D_list{2}))))<1);

id_hasPI3K = find(Mat_MODEL_mut(:,length(list_outputs)+2)==1);
id_hasRAS = find(Mat_MODEL_mut(:,length(list_outputs)+1)==1);

id_hasnoPI3K = find(Mat_MODEL_mut(:,length(list_outputs)+2)==0);
id_hasnoRAS = find(Mat_MODEL_mut(:,length(list_outputs)+1)==0);

% bdt_lowErbB2 = TreeBagger(500,Mat_MODEL_mut(id_lowErbB2,[1:19]),Mat_MODEL_mut(id_lowErbB2,end),'SampleWithReplacement','on','Method','classification','oobvarimp','on','MinLeaf',5);
% bdt_highErbB2 = TreeBagger(500,Mat_MODEL_mut(id_highErbB2,[1:19]),Mat_MODEL_mut(id_highErbB2,end),'SampleWithReplacement','on','Method','classification','oobvarimp','on','MinLeaf',5);

figure
scatter(log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasRAS)),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasRAS)),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasRAS)),find(ismember(list_outputs,D_list{5}))))+0.1)*8,zeros(length(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasRAS)),1)),1),'filled')
hold on
scatter(log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasRAS)),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasRAS)),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasRAS)),find(ismember(list_outputs,D_list{5}))))+0.1)*8,ones(length(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasRAS)),1)),1),'filled')
colormap(jet)
set(gcf,'Color','w')
xlabel('pS6 AUC')
ylabel('pAKT AUC')

figure
scatter(log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS)),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS)),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS)),find(ismember(list_outputs,D_list{5}))))+0.1)*8,zeros(length(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS)),1)),1),'filled')
hold on
scatter(log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS)),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS)),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS)),find(ismember(list_outputs,D_list{5}))))+0.1)*8,ones(length(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS)),1)),1),'filled')
colormap(jet)
set(gcf,'Color','w')
xlabel('pS6 AUC')
ylabel('pAKT AUC')

% figure
% scatter(log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS) & ismember(id_false,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS) & ismember(id_false,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS) & ismember(id_false,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{5}))))+0.1)*5,zeros(length(Mat_MODEL_mut(id_false(ismember(id_false,id_lowErbB2) & ismember(id_false,id_hasnoRAS) & ismember(id_false,find(Mat_MODEL_mut(:,20)==1))),1)),1),'filled')
% hold on
% scatter(log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS) & ismember(id_true,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{3})))),log(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS) & ismember(id_true,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{4})))),(1.1.^(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS) & ismember(id_true,find(Mat_MODEL_mut(:,20)==1))),find(ismember(list_outputs,D_list{5}))))+0.1)*5,ones(length(Mat_MODEL_mut(id_true(ismember(id_true,id_lowErbB2) & ismember(id_true,id_hasnoRAS) & ismember(id_true,find(Mat_MODEL_mut(:,20)==1))),1)),1),'filled')
% colormap(jet)
% set(gcf,'Color','w')
% xlabel('pS6 AUC')
% ylabel('pAKT AUC')

figure
scatter(log(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasPI3K)),find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasPI3K)),find(ismember(list_outputs,D_list{2})))),((Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasPI3K)),find(ismember(list_outputs,D_list{6}))))+1)*13,zeros(length(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasPI3K)),1)),1),'filled')
hold on
scatter(log(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasPI3K)),find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasPI3K)),find(ismember(list_outputs,D_list{2})))),((Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasPI3K)),find(ismember(list_outputs,D_list{6}))))+1)*13,ones(length(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasPI3K)),1)),1),'filled')
colormap(jet)
set(gcf,'Color','w')
xlabel('EGFR homodimerization')
ylabel('EGFR-ErbB2 heterodimerization')

figure
scatter(log(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasnoPI3K)),find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasnoPI3K)),find(ismember(list_outputs,D_list{2})))),((Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasnoPI3K)),find(ismember(list_outputs,D_list{6}))))+1.)*13,zeros(length(Mat_MODEL_mut(id_false(ismember(id_false,id_highErbB2) & ismember(id_false,id_hasnoPI3K)),1)),1),'filled')
hold on
scatter(log(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasnoPI3K)),find(ismember(list_outputs,D_list{1})))),log(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasnoPI3K)),find(ismember(list_outputs,D_list{2})))),((Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasnoPI3K)),find(ismember(list_outputs,D_list{6}))))+1.)*13,ones(length(Mat_MODEL_mut(id_true(ismember(id_true,id_highErbB2) & ismember(id_true,id_hasnoPI3K)),1)),1),'filled')
colormap(jet)
set(gcf,'Color','w')
xlabel('EGFR homodimerization')
ylabel('EGFR-ErbB2 heterodimerization')
