% can be used to merge plots of multiple data sets
% 
% usage (here jm is the model index):
% 
% arLoadData('data1', jm);
% arLoadData('data2', jm);
% ...
% arLoadData('dataN', jm);
% arMergePlotMulti(jm, {'data1', 'data2', ..., 'dataN'}, {'label for data1', 'label for data2', ..., 'label for dataN'}, 'new condition name');
% 
% comment: actually the names of the plots as defined in ar.model.plot.name
% are used

function arMergePlotMulti(jm, index_names, label_names, new_condition_name)

global ar

if(isempty(index_names))
    return
end

nNames = length(index_names);
index = nan(1,nNames);
conditions = cell(1,nNames);

for inames = 1:nNames

    if(~iscell(index_names))
        index(inames) = index_names(inames);
    else
        index_names{inames} = strrep(strrep(strrep(strrep(index_names{inames},'=','_'),'.',''),'-','_'),'/','_');
        index(inames) = -1;
        for jp = 1:length(ar.model(jm).plot)
            if(strcmp(ar.model(jm).plot(jp).name, index_names{inames}) && index(inames)==-1)
               index(inames) = jp;
            elseif(strcmp(ar.model(jm).plot(jp).name, index_names{inames}) && index(inames)~=-1)
                error('multiple matches for data set with name %s', index_names{inames});
            end
        end
        if(index(inames)==-1)
            error('could not find data set with name %s', index_names{inames});
        end
    end
    
end

% condition labels
for ic = 1:nNames
    
    conditions{ic} = cell(1,length(ar.model(jm).plot(index(ic)).dLink));
    for j=1:length(ar.model(jm).plot(index(ic)).dLink)
        if(isempty(ar.model(jm).plot(index(ic)).condition) || isempty(ar.model(jm).plot(index(ic)).condition{j}))
            conditions{ic}{j} = label_names{ic};
        else
            conditions{ic}{j} = [ar.model(jm).plot(index(ic)).condition{j} ' - ' label_names{ic}];
        end
    end
    
end

% response parameters
for ip = 2:nNames
    if(ar.model(jm).plot(index(1)).doseresponse==1 && ar.model(jm).plot(index(ip)).doseresponse==1)
        resppar1 = ar.model(jm).data(ar.model(jm).plot(index(1)).dLink(1)).response_parameter;
        resppar2 = ar.model(jm).data(ar.model(jm).plot(index(ip)).dLink(1)).response_parameter;
        if(~strcmp(resppar1, resppar2))
            ar.model(jm).plot(index(1)).response_parameter = 'dose';
            ar.model(jm).plot(index(ip)).response_parameter = 'dose';
        else
            ar.model(jm).plot(index(1)).response_parameter = resppar1;
            ar.model(jm).plot(index(ip)).response_parameter = resppar2;
        end
    elseif(ar.model(jm).plot(index(1)).doseresponse==1 || ar.model(jm).plot(index(ip)).doseresponse==1)
        error('trying to merge dose response with time course plots');
    end
end

% do merge
ar.model(jm).plot(index(1)).condition = conditions{1};
for ip = 2:nNames
    
    ar.model(jm).plot(index(1)).dLink = [ar.model(jm).plot(index(1)).dLink ar.model(jm).plot(index(ip)).dLink];
    ar.model(jm).plot(index(1)).condition = [ar.model(jm).plot(index(1)).condition conditions{ip}];
    if(exist('new_condition_name','var'))
        ar.model(jm).plot(index(1)).name = new_condition_name;
    else
        ar.model(jm).plot(index(1)).name = [ar.model(jm).plot(index(1)).name '_' ar.model(jm).plot(index(ip)).name];
    end
end

ar.model(jm).plot = ar.model(jm).plot(setdiff(1:length(ar.model(jm).plot),index(2:end)));
