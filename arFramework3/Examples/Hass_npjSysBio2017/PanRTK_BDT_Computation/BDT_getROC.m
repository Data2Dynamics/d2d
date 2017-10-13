%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   BDT_bootstrap_cluster                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% part of BDT_bootstrap script to obtain prediction Matrix and ROC curve
% for tree 'b' with testing data having features X and outcome Y
%
% ConfusionMatrix has architecture (0 = no growth, 1 = growth)
% ROWS=TRUE VALUES / COLS= PRED VALUES
%    0    1
% 0 TN   FP
% 1 FN   TP

function [xVal, yVal, auc, confMat, bdt_classes, predClass] = BDT_getROC(b,X_test,Y_test)

%get the predictions
[predClass, ClassScore] = b.predict(X_test);
if(iscell(predClass))
    [confMat,bdt_classes] = confusionmat(Y_test, str2double(predClass),'order',[0 1]);
    predClass = str2double(predClass);
else
    [confMat,bdt_classes] = confusionmat(Y_test, predClass,'order',[0 1]);   
end

%This part is currently not needed, enables to plot ROC curves

Y_names = unique(Y_test);
if(iscell(b.ClassNames))
    [Class_sort, I_class] = sort(str2double(b.ClassNames));    
    col_names = Class_sort;   
else
    [Class_sort, I_class] = sort(b.ClassNames);  
    col_names = Class_sort;
end
ClassScore_tmp = ClassScore;
for i = 1:size(ClassScore,2)
    ClassScore(:,i) = ClassScore_tmp(:,I_class(i));
end
hasClass = ismember(col_names,Y_names);

if(length(Y_names)<2)
    xVal=NaN(10,size(ClassScore,2));
    yVal=NaN(10,size(ClassScore,2));
    auc = NaN(1,size(ClassScore,2));
    confMat = NaN(size(ClassScore,2),size(ClassScore,2));
    return;
end
for j = 1:size(ClassScore,2)
    if(hasClass(j)==0)
        xVal_tmp = NaN;
        yVal_tmp = NaN;
        auc_tmp = NaN;
    else
        %get ROC curve and save it. Bit ugly, too, since size of vectors varies
        [xVal_tmp, yVal_tmp,~,auc_tmp] = perfcurve(Y_test, ClassScore(:,j),col_names(j));
    end
    if(j==1)
        xVal = xVal_tmp;
        yVal = yVal_tmp;
        auc = auc_tmp;    
    end
end