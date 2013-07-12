% can be used to merge plot of two data sets
% 
% usage (here jm=1 is the model index):
% 
% arLoadData('data1', jm);
% arLoadData('data2', jm);
% arMergePlot(jm, 'data1', 'data2', 'label for data1', 'label for data2');

function arMergePlot(jm, index_name1, index_name2, label_name1, label_name2)

global ar

if(isnumeric(index_name1))
    index1 = index_name1;
else
    index_name1 = strrep(strrep(strrep(strrep(index_name1,'=','_'),'.',''),'-','_'),'/','_');
    index1 = -1;
    for jp = 1:length(ar.model(jm).plot)
        if(strcmp(ar.model(jm).plot(jp).name, index_name1) && index1==-1)
           index1 = jp;
        elseif(strcmp(ar.model(jm).plot(jp).name, index_name1) && index1~=-1)
            error('multiple matches for data set with name %s', index_name1);
        end
    end
    if(index1==-1)
        error('could not find data set with name %s', index_name1);
    end
end

if(isnumeric(index_name2))
    index2 = index_name2;
else
    index_name2 = strrep(strrep(strrep(strrep(index_name2,'=','_'),'.',''),'-','_'),'/','_');
    index2 = -1;
    for jp = 1:length(ar.model(jm).plot)
        if(strcmp(ar.model(jm).plot(jp).name, index_name2) && index2==-1)
           index2 = jp;
        elseif(strcmp(ar.model(jm).plot(jp).name, index_name2) && index2~=-1)
            error('multiple matches for data set with name %s', index_name2);
        end
    end
    if(index2==-1)
        error('could not find data set with name %s', index_name2);
    end
end

ar.model(jm).plot(index1).dLink = [ar.model(jm).plot(index1).dLink ar.model(jm).plot(index2).dLink];
ar.model(jm).plot(index1).condition = {label_name1 label_name2};
ar.model(jm).plot(index1).name = [ar.model(jm).plot(index1).name '_' ar.model(jm).plot(index2).name];

ar.model(jm).plot = ar.model(jm).plot([1:(index2-1) (index2+1:length(ar.model(jm).plot))]);
