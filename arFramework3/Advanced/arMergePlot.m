% can be used to merge plot of two data sets
% 
% usage (here jm=1 is the model index):
% 
% arLoadData('data1', jm);
% arLoadData('data2', jm);
% arMergePlot(jm, 'data1', 'data2', 'label for data1', 'label for data2', 'new condition name');

function arMergePlot(jm, index_name1, index_name2, label_name1, label_name2, new_condition_name)

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

% condition labels
conditions1 = cell(1,length(ar.model(jm).plot(index1).dLink));
for j=1:length(ar.model(jm).plot(index1).dLink)
    if(isempty(ar.model(jm).plot(index1).condition) || isempty(ar.model(jm).plot(index1).condition{j}))
        conditions1{j} = label_name1;
    else
        conditions1{j} = [ar.model(jm).plot(index1).condition{j} ' & ' label_name1];
    end
end
conditions2 = cell(1,length(ar.model(jm).plot(index2).dLink));
for j=1:length(ar.model(jm).plot(index2).dLink)
    if(isempty(ar.model(jm).plot(index2).condition) || isempty(ar.model(jm).plot(index2).condition{j}))
        conditions2{j} = label_name2;
    else
        conditions2{j} = [ar.model(jm).plot(index2).condition{j} ' & ' label_name2];
    end
end

% response parameters
if(ar.model(jm).plot(index1).doseresponse==1 && ar.model(jm).plot(index2).doseresponse==1)
    resppar1 = ar.model(jm).data(ar.model(jm).plot(index1).dLink(1)).response_parameter;
    resppar2 = ar.model(jm).data(ar.model(jm).plot(index2).dLink(1)).response_parameter;
    if(~strcmp(resppar1, resppar2))
        ar.model(jm).plot(index1).response_parameter = 'dose';
        ar.model(jm).plot(index2).response_parameter = 'dose';
    else
        ar.model(jm).plot(index1).response_parameter = resppar1;
        ar.model(jm).plot(index2).response_parameter = resppar2;
    end
elseif(ar.model(jm).plot(index1).doseresponse==1 || ar.model(jm).plot(index2).doseresponse==1)
    error('trying to merge dose response with time course plots');
end

% do merge
ar.model(jm).plot(index1).dLink = [ar.model(jm).plot(index1).dLink ar.model(jm).plot(index2).dLink];
ar.model(jm).plot(index1).condition = [conditions1 conditions2];
if(exist('new_condition_name','var'))
    ar.model(jm).plot(index1).name = new_condition_name;
else
    ar.model(jm).plot(index1).name = [ar.model(jm).plot(index1).name '_' ar.model(jm).plot(index2).name];
end 
ar.model(jm).plot = ar.model(jm).plot([1:(index2-1) (index2+1:length(ar.model(jm).plot))]);
