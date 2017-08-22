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

function [xVal, yVal, auc, confMat, bdt_classes, predClass] = BDT_getROC(b,X_test,Y_test)

%get the predictions
[predClass, ClassScore] = b.predict(X_test);
if(iscell(predClass))
    confMat = confusionmat(Y_test, str2double(predClass));
    bdt_classes = unique(str2double(predClass));
    predClass = str2double(predClass);
else
    confMat = confusionmat(Y_test, predClass);
    bdt_classes = unique(predClass);

end
Y_names = unique(Y_test);%[-1 0 1];
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
    else
        if(length(xVal_tmp)~=size(xVal,1))
            if(length(xVal_tmp)<size(xVal,1))
               xVal_tmp = [xVal_tmp;NaN(size(xVal,1)-length(xVal_tmp),1)];
               yVal_tmp = [yVal_tmp;NaN(size(xVal,1)-length(yVal_tmp),1)];
            else
               xVal = [xVal;NaN(length(xVal_tmp)-size(xVal,1),size(xVal,2))];
               yVal = [yVal;NaN(length(yVal_tmp)-size(yVal,1),size(yVal,2))];
            end
        end
        xVal(:,j)=xVal_tmp;
        yVal(:,j)=yVal_tmp;
        auc(j) = auc_tmp;
    end
end