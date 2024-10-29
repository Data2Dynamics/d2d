% arSetTRange(model, data, [Tvar], [minT], [maxT])
%
% Set of the measurement times T in model definition file. Used for RTFs.
%
%   model   path to model definition file or model name
%   data    path to data file or data name
%   Tvar    name of variable in model definition file
%           ['T']
%   minT    minimum time for time range T
%           'data' = load data and find min(T)
%           [0]
%   maxT    maximum time for time range T
%           'data' = load data and find max(T)
%           [max(t in data)]
%

function T = arSetTRange(model, data, Tvar, minT, maxT)
    
    if isempty(strmatch(strtok(model, '/'),'Models','exact')) % check if model input is model name or path
        model = ['Models/',model];
    end
    modelSplitDot = split(model, '.');
    if isempty(strmatch(modelSplitDot{end},'def','exact')) 
        model = [model,'.def'];
    end
    
    if isempty(strmatch(strtok(data, '/'),'Data','exact')) % check if data input is data name or path
        data = ['Data/',data];
    end
    dataSplitDot = split(data, '.');
    if isempty(strmatch(dataSplitDot{end},'xlsx','exact'))
        data = [data,'.xlsx'];
    end

    [d,s] = xlsread(data); % read data file
    index_time = find(contains(s,'time'));

    if ~exist("Tvar","var") || isempty(Tvar) % check if Tvar is given
        Tvar = 'T';
    end
    
    if ~exist("minT","var") || isempty(minT) % check if minT is given
        minT = 0;
    elseif strcmp(minT,'data')    
        minT = min(d(:,index_time));
    elseif isnumeric(minT)
        minT = minT;
    else
        error('minT must be numeric or ''data''')
    end

    if ~exist("maxT","var") || isempty(maxT) % check if maxT is given
        maxT = max(d(:,index_time));
    elseif strcmp(maxT,'data')    
        maxT = max(d(:,index_time));
    elseif isnumeric(maxT)
        maxT = maxT;
    else
        error('maxT must be numeric or ''data''')
    end
    
    % calculate T as range from minT to maxT
    T = maxT-minT;
    
    % read model def file and replace T with new value
    Tvar = 'T';
    filetext  = fileread(model);
    exp_search = ['\n',Tvar,'[^\n]*'];
    match = regexp(filetext,exp_search,'match');
    expr = match{1};
    expr_new = ['\n',Tvar,'   "',num2str(T),'"'];
    filetext_new = strrep(filetext,expr,expr_new);
    
    % write new model def file
    fid = fopen(model,'w');
    fprintf(fid, filetext_new);
    fclose(fid);
end