%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   Setup                                   %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%adapted from arLoadData of arFramework (www.data2dynamics.org)
% Loads data file into the first cell line in data struct and fills up the
% respective conditions, sets ar.fit to 0 for all others and creates new ar.model(1).plot struct.

function LoadABData(name, extension)

global ar

% load model from mat-file
if(~exist('Data','dir'))
    error('folder Data/ does not exist')
end
if(~exist(['Data/' name '.def'],'file'))
    if(~exist(['Data/' name '.xls'],'file') && ~exist(['Data/' name '.csv'],'file') && ~exist(['Data/' name '.xlsx'],'file'))
        error('data definition file %s.def does not exist in folder Data/', name)
    else
%         arFprintf(1, '\ncreating generic .def file for Data/%s ...\n', name);
%         copyfile(which('data_template.def'),['./Data/' name '.def']);
        error('create def file!!')
    end
else
    if(~exist(['Data/' name '.xls'],'file') && ~exist(['Data/' name '.csv'],'file') && ~exist(['Data/' name '.xlsx'],'file'))
        warning('data file corresponding to %s.def does not exist in folder Data/', name)
    end
end
    
if(~exist('m','var') || isempty(m))
    m = length(ar.model);
end
if(exist('m','var') && ischar(m))
    for jm=1:length(ar.model)
        if(strcmp(m, ar.model(jm).name))
            m = jm;
        end
    end
    if(ischar(m))
        error('Model %s was not found', m);
    end
end

if(exist('extension','var') && isnumeric(extension) && ...
        ~(isempty(extension) && nargin>3))
    error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
        'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
end


if(~exist('extension','var') || isempty(extension))
    extension = 'xls';
end
if(~exist('removeEmptyObs','var'))
    removeEmptyObs = false;
else
    if(ischar(removeEmptyObs))
        error(['arLoadData(name, m, d, ...) input argument d is deprecated !!! ' ...
            'Please see new usage arLoadModel(name, m, extension, removeEmptyObs) and function help text.']);
    end
end

% XLS file
if(~strcmp(extension,'none') && ( ...
    (exist(['Data/' name '.xlsx'],'file') && strcmp(extension,'xls')) ||...
    (exist(['Data/' name '.xls'],'file') && strcmp(extension,'xls')) || ...
    (exist(['Data/' name '.csv'],'file') && strcmp(extension,'csv'))))
    arFprintf(2, 'loading data from file Data/%s.%s...\n', name, extension);

    % read from file
    if(strcmp(extension,'xls'))
        warntmp = warning;
        warning('off','all')
        
        if (exist(['Data/' name '.xls'],'file'))      
            [data, Cstr] = xlsread(['Data/' name '.xls']);
        elseif (exist(['Data/' name '.xlsx'],'file'))      
            [data, Cstr] = xlsread(['Data/' name '.xlsx']);
        end
        
        if(length(data(1,:))>length(Cstr(1,:)))
            data = data(:,1:length(Cstr(1,:)));
        end
        
        warning(warntmp);
        
        header = Cstr(1,2:end);
        header = strrep(header,' ',''); % remove spaces which are sometimes in the column header by accident    
        times = data(:,1);
        qtimesnonnan = ~isnan(times);
        times = times(qtimesnonnan);
        data = data(qtimesnonnan,2:end);
        if(size(data,2)<length(header))
            data = [data nan(size(data,1),length(header)-size(data,2))];
        end
        
        Cstr = Cstr(2:end,2:end);
        dataCell = cell(size(data));
        for j1 = 1:size(data,1)
            for j2 = 1:size(data,2)
                if(isnan(data(j1,j2)))
                    if(j1<=size(Cstr,1) && j2<=size(Cstr,2) && ~isempty(Cstr{j1,j2}))
                        dataCell{j1,j2} = Cstr{j1,j2};
                    else
                        dataCell{j1,j2} = header{j2};
                    end
                else
                    dataCell{j1,j2} = num2str(data(j1,j2));
                end
            end
        end
        
    elseif(strcmp(extension,'csv'))
        [header, data, dataCell] = arReadCSVHeaderFile(['Data/' name '.csv'], ',', true);

        header = header(2:end);
        times = data(:,1);
        data = data(:,2:end);
        dataCell = dataCell(:,2:end);

    end
    nd = length(ar.model(1).plot(1).dLink);
    conditions = NaN(nd,1);
    cond_names = cell({});
    for j=1:nd
        for i=1:length(ar.model(1).data(ar.model(1).plot(1).dLink(j)).condition)
            if(strfind(ar.model(1).data(ar.model(1).plot(1).dLink(j)).condition(i).parameter,'_level')>0)
                conditions(j,i)=str2double(ar.model(1).data(ar.model(1).plot(1).dLink(j)).condition(i).value);
                if(j==1)
                    cond_names = [cond_names ar.model(1).data(ar.model(1).plot(1).dLink(j)).condition(i).parameter];
                end
            end
        end   
    end 
    ar.model(m).plot(end+1).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
    ar.model(m).plot(end).doseresponse = false;
    ar.model(m).plot(end).doseresponselog10xaxis = true;
    ar.model(m).plot(end).ny = length(ar.model(m).data(1).y);
    ar.model(m).plot(end).condition = {};
    ar.model(m).plot(end).dLink = [];
    ar.model(m).qPlotYs = [ar.model(m).qPlotYs 1];
    ar.model(m).qPlotVs = [ar.model(m).qPlotVs 0];
    ar.model(m).qPlotXs = [ar.model(m).qPlotXs 0];
    
    jplot = length(ar.model(m).plot);

    % collect parameters conditions
    pcond = union(ar.model(m).data(1).p, ar.model(m).data(1).pcond); %R2013a compatible
    [ar, d] = setConditions(ar, m, name, header, times, data, dataCell, pcond, conditions, cond_names, jplot);

    % Check whether the user specified any variables with reserved words.
    checkReserved(m, d);
    
else
    warning('Cannot find data file corresponding to %s', name);
    
end

ar = orderfields(ar);
ar.model = orderfields(ar.model);
ar.model(m).data = orderfields(ar.model(m).data);
ar.model(m).plot = orderfields(ar.model(m).plot);

function checkReserved(m, d)
    global ar;

    % Check whether the user specified any variables with reserved words.
    for a = 1 : length( ar.model(m).data(d).fu )
        arCheckReservedWords( symvar(ar.model(m).data(d).fu{a}), sprintf( 'input function of %s', ar.model(m).data(d).name ), ar.model(m).u{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fy )
        arCheckReservedWords( symvar(ar.model(m).data(d).fy{a}), sprintf( 'observation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
    for a = 1 : length( ar.model(m).data(d).fystd )
        arCheckReservedWords( symvar(ar.model(m).data(d).fystd{a}), sprintf( 'observation standard deviation function of %s', ar.model(m).data(d).name ), ar.model(m).data(d).y{a} );
    end
%     for a = 1 : length( ar.model(m).data(d).fp )
%         arCheckReservedWords( symvar(ar.model(m).data(d).fp{a}), sprintf( 'condition parameter transformations of %s', ar.model(m).data(d).name ), ar.model(m).data(d).p{a} );
%     end   
    arCheckReservedWords( ar.model(m).data(d).p, 'parameters' );
    arCheckReservedWords( ar.model(m).data(d).y, 'observable names' );

function [ar,d] = setConditions(ar, m, name, header, times, data, dataCell, pcond, conditions, cond_names, jplot)

% matVer = ver('MATLAB');

% normalization of columns
nfactor = max(data, [], 1);

qobs = ismember(header, ar.model(m).data(1).y) & sum(~isnan(data),1)>0; %R2013a compatible
qhasdata = ismember(ar.model(m).data(1).y, header(qobs)); %R2013a compatible

% conditions
qcond = ismember(header, pcond); %R2013a compatible

condi_header = header(qcond);
q_level = ~cellfun(@isempty,(strfind(condi_header,'_level')));
q_level = find(q_level==1);

if ~isempty(dataCell)
    [condis, ~, jcondis] = uniqueRowsCA(dataCell(:,qcond));
else
    [condis, ~, jcondis] = unique(data(:,qcond),'rows');
    condis = mymat2cell(condis);
end

active_condi = false(size(condis(1,:)));
tmpcondi = condis(1,:);
for j1=2:size(condis,1)
    for j2=1:size(condis,2)
        active_condi(j2) = active_condi(j2) | (~strcmp(tmpcondi{j2}, condis{j1,j2}));
    end
end

for j=1:size(condis,1)

    arFprintf(2, 'local condition #%i:\n', j)        
    
    cond_tmp = NaN(size(condis(j,length(q_level))));
    for jj=1:length(q_level)
        cond_tmp(jj) = str2double(condis{j,q_level(jj)});        
    end
    d = strmatch(cond_tmp,conditions);
    if(isempty(d))
        continue;
    end
    % initial setup
    ar.model(m).data(d).name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
    ar.model(m).data(d).uNames = {};    
    for jj=1:size(condis,2)
    % plot
        if(active_condi(jj))
            if(ar.model(m).data(d).doseresponse==0 || ~strcmp(condi_header{jj}, ar.model(m).data(d).response_parameter))
                if(length(ar.model(m).plot(jplot).condition) >= j && ~isempty(ar.model(m).plot(jplot).condition{j}))
                    ar.model(m).plot(jplot).condition{j} = [ar.model(m).plot(jplot).condition{j} ' & ' ...
                        ar.model(m).data(d).condition(jj).parameter '=' ...
                        ar.model(m).data(d).condition(jj).value];
                else
                    ar.model(m).plot(jplot).condition{j} = [ar.model(m).data(d).condition(jj).parameter '=' ...
                        ar.model(m).data(d).condition(jj).value];
                end
            end
        end
    end
    
    ar.model(m).plot(jplot).dLink(end+1) = d;
    
    qvals = jcondis == j;
    ar = setValues(ar, m, d, header, nfactor, data(qvals,:), times(qvals));
    ar.model(m).data(d).tLim(2) = round(max(times)*1.1);             

end


function C = mymat2cell(D)
C = cell(size(D));
for j=1:size(D,1)
    for jj=1:size(D,2)
        C{j,jj} = num2str(D(j,jj));
    end
end



function ar = setValues(ar, m, d, header, nfactor, data, times)
ar.model(m).data(d).tExp = times;
ar.model(m).data(d).yExp = nan(length(times), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = nan(length(times), length(ar.model(m).data(d).y));

for j=1:length(ar.model(m).data(d).y)
    q = ismember(header, ar.model(m).data(d).y{j}); %R2013a compatible
    
    if(sum(q)==1)
        ar.model(m).data(d).yExp(:,j) = data(:,q);
        arFprintf(2, '\t%20s -> %4i data-points assigned', ar.model(m).data(d).y{j}, sum(~isnan(data(:,q))));
        
        % normalize data
        if(ar.model(m).data(d).normalize(j))
            ar.model(m).data(d).yExp(:,j) = ar.model(m).data(d).yExp(:,j) / nfactor(q);
            arFprintf(2, ' normalized');
        end
        
        % log-fitting
        if(ar.model(m).data(d).logfitting(j))
            qdatapos = ar.model(m).data(d).yExp(:,j)>0;
            ar.model(m).data(d).yExp(qdatapos,j) = log10(ar.model(m).data(d).yExp(qdatapos,j));
            ar.model(m).data(d).yExp(~qdatapos,j) = nan;
            if(sum(~qdatapos)==0)
                arFprintf(2, ' for log-fitting');
            else
                arFprintf(2, ' for log-fitting (%i values <=0 removed)', sum(~qdatapos));
            end
        end
        
        % empirical stds
        qstd = ismember(header, [ar.model(m).data(d).y{j} '_std']); %R2013a compatible
        if(sum(qstd)==1)
            ar.model(m).data(d).yExpStd(:,j) = data(:,qstd);
            arFprintf(2, ' with stds');
            if(ar.model(m).data(d).normalize(j))
                ar.model(m).data(d).yExpStd(:,j) = ar.model(m).data(d).yExpStd(:,j) / nfactor(q);
                arFprintf(2, ' normalized');
            end
        elseif(sum(qstd)>1)
            error('multiple std colums for observable %s', ar.model(m).data(d).y{j})
        end
        
    elseif(sum(q)==0)
        arFprintf(2, '*\t%20s -> not assigned', ar.model(m).data(d).y{j});
    else
        error('multiple data colums for observable %s', ar.model(m).data(d).y{j})
    end
    
    arFprintf(1, '\n');
end

