% D = arCreateDataStructure(fieldsvar, mixingCond, [name], [opts])
%
%
% *** ATTENTION ***
% Use this function cautiously! It could overrides the existing data struct.
%
% This function can be used to create data structs before compiling. This
% data struct can be added to the existing data struct (e.g. the one read
% from def file) it or it can be used as an independent data struct. for
% more information about using this function refer to the arLoadData
% 
%
%   fieldsvar            cell     information provided by user
%   mixingCond           double   mixing condition = [x x x]
%                                 elemenst are either zero or on
%                                 1st element: data struct
%                                       0: discars new data stuct
%                                          othe elements are not important
%                                       1: takes into account the data struct  
%                                 2nd element: def information
%                                       0: updates the data struct by the
%                                          new one
%                                       1: replace the existing data struct
%                                          with the new one 
%                                 3rd element: data file
%                                       0: updates the data file by the
%                                          new one
%                                       1: replace the existing data file
%                                          with the new one
%   name                 char     if prefer to change the default name
%                                 in data struct
%                                 ['arEmptyDataStruct']
%   opts                 strucr   to extract data path
%                                 ['Data/']
%
%   available fields in .def file
%   DESCRIPTION: description
%   PREDICTOR: response_parameter, t, tUnits, tLim, tExpLim
%   INPUTS: u, fu, uNames
%   OBSERVABLES: y, yUnits, logfitting, logplotting, normalize, fy, yNames
%   ERRORS: y, fystd
%   CONDITIONS: pcond, fpcond
%   RANDOM: prand, rand_type
%   data: (tExp, yExp, yExpStd) or (header, data)
%   (SUBSTITUTIONS & PARAMETERS will be added)
%
%   Example:
%   description = {'user input data struct' ; '2nd line'};
%   tLim = [0 10];
%   y = {'y1' 'y2' 'y3'};
%   fy = {'fy1'; 'fy2'; 'fy3'};
%   fystd = {'fystd1'; 'fystd2'; 'fystdy3'};
%   yUnits =  {'C' 'au' 'conc'}
%   yUnits = repmat(yUnits,[3,1]);
%   tExp=[1 2 3 4];
%   yExp=[11 12 13; 21 22 23; 31 32 33; 41 42 43 44];
%   yExp=[11 12 13; 21 22 23; 31 32 33; 41 42 43 44];
%   pcond = {'cond1'};
%   fpcond = {'fpcond'};
%   fielsvar = {'description' description 'tLim' tLim 'y' y 'fy' fy 'fystd' fystd' 'yUnits' yUnits 'pcond' pcond 'fpcond' fpcond};
%
%   If you have a dose-response parameter you can use the following method:
%   you should define y, data and dose-response information exactly like xls, xlsx or cvs file
%   header = {'time' 'dose-response' 'y1' 'y2' 'y3' 'y1_std' 'y2_std' 'y3_std'};
%   data = nx8 matrix
%   you may use this method even for non dose-response problems as well.
 

function D = arCreateDataStructure(fieldsvar, mixingCond, name, opts)


if ~exist('name','var') || ~ischar(name)
    name='';
end
name = strrep(name,' ','');

if ~exist('opts','var') || isempty(opts) || isempty(opts.datapath_args)
    DataPath = 'Data/';
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end

if ~exist('fieldsvar','var') || isempty(fieldsvar)
    D = arEmptyDataStruct;
    D.DataPath = DataPath;
    D = checkDef(name, D, true);
    return
else
    fprintf('\nParsing data struct introduced by user...');
end


isDef = true;   % data struct should contains all the essentials def file infomation
%isData = true;  % experimental data should be self-consistent
if exist('mixingCond','var') && ~isempty(mixingCond) && isnumeric(mixingCond)
    if mixingCond(2)==0
        isDef = false;
    end
else
    error('Invalid mixingCond')
end


D = arEmptyDataStruct;

D.header = cell(0);
D.data = [];
D.xlstimevar = cell(0);
D.xlstimes = [];
D.xlsheader = cell(0);
D.xlsdata = [];
D.xlsdataCell = cell(0);

D.DataPath = DataPath;

D = extractData(D,fieldsvar);

D = checkDef(name, D, isDef);

D = checkData(D);

D.tExp=[];
D.yExp=[];
D.yExpStd=[];
D=rmfield(D,'data');
D=rmfield(D,'header');

D = orderfields(D);

fprintf('[OK]\n');


end


function D = checkData(D)
if isempty(D.tLimExp)
    D.tLimExp = [nan nan];
end
    
%Performing checks (if experimental data is provided through fieldsvar)
if ~isempty(D.tExp)
    if sum(isnan(D.tExp))>0
        error('NaN in tExp')
    elseif ~any(size(D.tExp)==1)
        error('D.tExt should be one-dimensional array');
    elseif size(D.yExp,2) ~= size(D.y,2)
        error('size(D.yExp,2) ~= size(D.y,2)');
    elseif size(D.yExp,1) ~= length(D.tExp)
        error('size(D.yExp,2) ~= size(D.y,2)');
    elseif any(size(D.yExpStd) - size(D.yExp))
        error('size(D.yExp) ~= size(D.yExpStd)');
    elseif length(D.tLimExp) ~= 2
        error('length(D.tLimExp)~=2');
    end
end
end

function D = checkDef(name,D,isDef)

% preliminary self-consistency check
if ~isempty(D.u) && (isempty(D.fu) || length(D.fu)~=length(D.u))
    error('inappropriate INPUTS definition');
end
if ~isempty(D.y) && (isempty(D.fy) || length(D.fy)~=length(D.y))
    error('inappropriate OBSERVABLES definition');
elseif ~isempty(D.y) && (isempty(D.fystd) || length(D.fystd)~=length(D.y))
    error('inappropriate OBSERVABLES definition');
end
if ~isempty(D.pcond) && (isempty(D.fpcond) || length(D.fpcond)~=length(D.pcond))
    error('inappropriate CONDITIONS definition');
end
if ~isempty(D.prand) && (isempty(D.rand_type) || length(D.rand_type)~=length(D.prand))
    error('inappropriate RANDOM definition');
end
if ~isempty(D.fu) && isempty(D.u)
    error('inappropriate INPUTS definition, you forgot to define D.u');
end
if ~isempty(D.fy) && isempty(D.y)
    error('inappropriate OBSERVABLES definition, you forgot to define D.y');
end
if ~isempty(D.fystd) && isempty(D.y)
    error('inappropriate OBSERVABLES definition, you forgot to define D.y');
end
if ~isempty(D.fpcond) && isempty(D.pcond)
    error('inappropriate CONDITIONS definition, you forgot to define D.pcond');
end
if ~isempty(D.rand_type) && isempty(D.prand)
    error('inappropriate RANDOM definition, you forgot to define D.prand');
end

if ~isempty(D.y)
    if isempty(D.yNames)
        D.yNames = D.y;
    elseif length(D.yNames)~=length(D.y)
        if length(D.yNames)<length(D.y)
            D.yNames(end+1:length(D.y))= D.y(length(D.yNames)+1:end);
        else
            error('length D.yNames > D.y')
        end
    end
end


if ~isempty(D.u)
    if isempty(D.uNames)
        D.uNames = D.u;
    elseif length(D.uNames)~=length(D.u)
        if length(D.uNames)<length(D.u)
            D.uNames(end+1:length(D.u))= D.u(length(D.uNames)+1:end);
        else
            error('length D.uNames > D.u')
        end
    end
end



if isDef
    
    if isempty(D.name)
        if ~isempty(name)
            D.name=[name '_AutoGen'] ;
        else
            D.name  = 'arEmptyDataStruct';
        end
    end
    
    
    if isfield(D,'DataPath')
        D.path=[pwd '/' D.DataPath];
        D = rmfield(D,'DataPath');
    end
    
    if isempty(D.description)
        D.description  = {'Created by arEmptyDataStruct'};
    end
    
    if isempty(D.doseresponse)
        D.doseresponse=0;     
    end
    
    % time
    if isempty(D.t)
        D.t='t';
    end
    if isempty(D.tUnits) % required for plotting
        D.tUnits = {'T'  'min'  'time'};
    end
    if isempty(D.tLim)
        if ~isempty(D.tExp)
            D.tLim = [min(D.tExp),max(D.tExp)];
        else
            D.tLim = [0 100];
        end
    end
    if sum(isnan(D.tLim))>0
        error('NaN in tLim');
    end
    if isempty(D.tLimExp)
        D.tLimExp = [NaN NaN]; % will be updated in arAddDataDefInfo
    end
    
    % observables
    if ~isempty(D.y)
        if isempty(D.yUnits) % required for plotting
            D.yUnits = repmat({'C'  'au'  'conc.'},size(D.y'));
        end
        if isempty(D.logplotting)
            D.logplotting = zeros(size(D.y));
        end
        if isempty(D.logfitting)
            D.logfitting = zeros(size(D.y));
        end
        if isempty(D.normalize)
            D.normalize = zeros(size(D.y));
        end
        
        if any(size(D.y)-size(D.logfitting))
            error('size(D.logfitting)~=size(D.y)');
        elseif any(size(D.y)-size(D.logplotting))
            error('size(D.logplotting)~=size(D.y)');
        elseif any(size(D.y)-size(D.normalize))
            error('size(D.normalize)~=size(D.y)');
        elseif any(size(D.y)-size(D.yNames))
            error('size(D.yNames)~=size(D.y)');
        elseif any(length(D.y)-size(D.yUnits,1))
            error('size(D.yUnits,1)~=size(D.y)');
        elseif size(D.yUnits,2)~=3
            error('size(D.yUnits,2)~=3');
        elseif size(D.yUnits,1)~=size(D.y,2)
            error('size(D.yUnits,1)~=size(D.y)');
        elseif any(length(D.y)-length(D.fystd))
            error('size(D.fystd)~=size(D.y)');
        end
        
    end
    
    
else
    
    if ~isempty(name)
        D.name=[name '_AutoGen'] ;
    end
    if ~isempty(D.tLim) && sum(isnan(D.tLim))>0
        error('NaN in tLim');
    end
    
    if ~isempty(D.y)        
        if ~isempty(D.logplotting) && length(D.logfitting)~=length(D.y)
            error('size(D.logfitting)~=size(D.y)');
        end
        if ~isempty(D.logfitting) && length(D.logfitting)~=length(D.y)
            error('size(D.logfitting)~=size(D.y)');
        end
        if ~isempty(D.normalize) && length(D.normalize)~=length(D.y)
            error('size(D.normalize)~=size(D.y)');
        end
        if ~isempty(D.yNames) && length(D.yNames)~=length(D.y)
            error('length(D.yNames)~=length(D.y)');
        end
        if ~isempty(D.yUnits) && size(D.yUnits,1)~=length(D.y)
            error('size(D.yUnits,1)~=length(D.y)');
        end
        if ~isempty(D.yUnits) && size(D.yUnits,2)~=3
            error('size(D.yUnits,2)~=3');
        end       
    end
         
end

end

function D = extractData(D, fieldsvar)


fn = fieldsvar(1:2:(end-1));
fv = fieldsvar(2:2:end);

if length(fn)~=length(fv)
    error('number of variables does not match number of values')
end

for i=1:length(fn)
    if isfield(D,fn{i})
        D.(fn{i})=fv{i};      
    else
        error('unknown field %s',fn{i})
    end       
end


if ~isempty(D.doseresponse) && ~isnumeric(D.doseresponse)
    error('doseresponse should be numeric')
end

% if ~isempty(D.response_parameter)
%     D.doseresponse = 1;
% end
   
if ~isempty(D.logfitting)
    if ~isnumeric(D.logfitting)
        error('logfitting should be numeric')
    elseif ~isrow(D.logfitting)
        D.logfitting=transpose(D.logfitting);
    end
end

if ~isempty(D.logplotting)
    if ~isnumeric(D.logplotting)
        error('logplotting should be numeric')
    elseif ~isrow(D.logplotting)
        D.logplotting=transpose(D.logplotting);
    end
end

if ~isempty(D.logplotting) && isempty(D.logfitting)
    D.logfitting=D.logplotting;
elseif isempty(D.logplotting) && ~isempty(D.logfitting)
    D.logplotting=D.logfitting;
end

if ~isempty(D.normalize)
    if ~isnumeric(D.normalize)
        error('normalize should be numeric')
    elseif ~isrow(D.normalize)
        D.normalize=transpose(D.normalize);
    end
end

if ~isempty(D.rand_type) && ~isnumeric(D.rand_type)
    error('rand_type should be numeric')
end

if ~isempty(D.pcond)
    if isrow(D.pcond)
       D.pcond=transpose(D.pcond); 
    end
end

if ~isempty(D.fpcond)
    if isrow(D.fpcond)
       D.fpcond=transpose(D.fpcond);
    end
    for i=1:length(D.fpcond)
        D.fpcond{i} = ['(' D.fpcond{i} ')'];
    end
end



% time
if ~isempty(D.t) && ~ischar(D.t)
    error('predictor should be char')
end

if ~isempty(D.tUnits) && length(D.tUnits)~=3
    error('length(D.tUnits)~=3')
    if ~isrow(D.tUnits)
        D.tUnits=transpose(D.tUnits);
    end
end

if ~isempty(D.tLim)
    if ~isnumeric(D.tLim)
        error('tLim should be numeric')
    elseif length(D.tLim)~=2
        error('length(D.tLim)~=2')
    elseif ~isrow(D.tLim)
        D.tLim=transpose(D.tLim);
    end
end

if ~isempty(D.tLimExp)
    if ~isnumeric(D.tLimExp)
        error('tLimExp should be numeric')
    elseif length(D.tLimExp)~=2
        error('length(D.tLimExp)~=2')
    elseif ~isrow(D.tLimExp)
        D.tLimExp=transpose(D.tLimExp);
    end
end

% inputs
if ~isempty(D.u)
    D.u = strrep(strrep(strtrim(D.u),'  ',''),' ','');
    if ~isrow(D.u)
        D.u=transpose(D.u);
    end
end
if ~isempty(D.fu)
    D.fu = strtrim(D.fu);
    if isrow(D.fu)
        D.fu=transpose(D.fu);
    end
end
if ~isempty(D.uNames)
    D.uNames = strtrim(D.uNames);
    if ~isrow(D.uNames)
        D.uNames=transpose(D.uNames);
    end
end

if ~isempty(D.yNames)
    D.yNames = strtrim(D.yNames);
    if ~isrow(D.yNames)
        D.yNames=transpose(D.yNames);
    end
end

if ~isempty(D.yUnits)
    D.yUnits = strtrim(D.yUnits);
end

if ~isempty(D.y)
    D.y = strrep(strrep(strtrim(D.y),'  ',''),' ','');
    if ~isrow(D.y)
        D.y=transpose(D.y);
    end
end

if ~isempty(D.fy)
    D.fy = strtrim(D.fy);
    if isrow(D.fy)
        D.fy=transpose(D.fy);
    end
end

if ~isempty(D.fystd)
    D.fystd = strtrim(D.fystd);
    if isrow(D.fystd)
        D.fystd=transpose(D.fystd);
    end
end

    

if isfield(D,'data') && ~isempty(D.data)
    
    if length(D.header)~=size(D.data,2)
        error('header length should be equal to the number of the columns in data');
    end
    
    timevar = D.header(1);
    times = D.data(:,1);
    header = D.header(2:end);
    data = D.data(:,2:end);
    %dataCell = cellfun(@num2str,num2cell(data),'UniformOutput',false);  
    
    if D.doseresponse==1
        if ismember(D.response_parameter, header)
            idx = ismember(header, D.response_parameter);
            if sum(isnan(data(:,idx)))>0
                error('does-response parameter should be numeric');
            end
            header_tmp = header(~idx);
            data_tmp = data(:,~idx);
        else
            error('dose-response parameter not found in data array');
        end
    else
        header_tmp = header;
        data_tmp = data;
    end
    
    % if rem(size(data,2),2)~=0
    %     error('number of yExp is not equal to number of yExpStd in data file');
    % end
    
    idx_y = ~endsWith(header_tmp,'_std');
    idx_yErr = endsWith(header_tmp,'_std');
    
    y = strrep(strtrim(header_tmp(idx_y)),' ','');
    y_err = strrep(strtrim(strrep(header_tmp(idx_yErr),'_std','')),' ','');
    
    if isempty(D.y)
        error('please input observables y');
    end
    if  ~all(ismember(y,D.y)==1) && ~all(ismember(D.y,y)==1)
        error('observables in header is not the same as y');
    end
      
    if length(y)~= length(y_err) || ~all(ismember(y,y_err)==1)
        error('number of yExp is not equal to the number of yExpStd in data file');
    end
    
    yExp = data_tmp(:,idx_y);
    yExpStd = data_tmp(:,idx_yErr);
    
    [q,qc]=ismember(y,y_err);
    yExpStd=yExpStd(:,qc);
        
    dataCell = cell(size(data));
    for i = 1:size(data,1)
        for j = 1:size(data,2)
            if isnan(data(i,j))
                dataCell{i,j} = header{j};
            else
                dataCell{i,j} = num2str(data(i,j));
            end
        end
    end    
    
    D.tExp = times;
    D.y = y;
    D.yExp = yExp;
    D.yExpStd = yExpStd;
    
    D.xlstimevar = timevar;
    D.xlstimes = times;
    D.xlsheader = header;
    D.xlsdata = data;
    D.xlsdataCell = dataCell;
    
    D.extension = 'userDefined';
    
    

elseif ~isempty(D.yExp) && isequal(D.doseresponse,1)    
    error('please use arrays data and header in case of dose-response, instead of tExp, yExp and yExpStd')    
else
    
    %creating raw data, like reading XLS!
    if ~isempty(D.tExp)
        if ~isnumeric(D.tExp)
            error('tExp should be numeric')
        else
            if isrow(D.tExp)
                D.tExp=transpose(D.tExp);
            end
        end
        D.extension = 'userDefined';
        D.xlstimes = D.tExp;
        D.xlstimevar = {'time'};
    end
    if ~isempty(D.yExp)
        D.xlsheader = D.y;
        if isnumeric(D.yExp)
            D.xlsdata=D.yExp;
            D.xlsdataCell = split(strtrim(cellstr(num2str(D.yExp))));
            for i = 1:size(D.xlsdataCell,1)
                for j = 1:size(D.xlsdataCell,2)
                    if strcmp(D.xlsdataCell{i,j},'NaN')
                        D.xlsdataCell(i,j) = D.xlsheader(j);
                    end
                end
            end
        end
    end
    if ~isempty(D.yExpStd)
        for i=1:length(D.y)
            tmp.xlsheader{i}=[D.y{i} ,'_std'];
        end
        if isnumeric(D.yExpStd)
            tmp.xlsdata=D.yExpStd;
            tmp.xlsdataCell = split(strtrim(cellstr(num2str(D.yExpStd))));
            for i = 1:size(tmp.xlsdataCell,1)
                for j = 1:size(tmp.xlsdataCell,2)
                    if strcmp(tmp.xlsdataCell{i,j},'NaN')
                        tmp.xlsdataCell(i,j) = tmp.xlsheader(j);
                    end
                end
            end
        end
        D.xlsheader = cat(2, D.xlsheader, tmp.xlsheader);
        D.xlsdata = cat(2, D.xlsdata, tmp.xlsdata);
        D.xlsdataCell = cat(2, D.xlsdataCell, tmp.xlsdataCell);
    end
 
end


end

function D = arEmptyDataStruct
    
                     % D = struct;
             D.condition = [];
           D.description = '';
          D.doseresponse = [];
                    D.fp = cell(0);
                D.fpcond = cell(0);
                    D.fu = cell(0);
                    D.fy = cell(0);
                 D.fystd = cell(0);
            D.logfitting = [];
           D.logplotting = [];
                  D.name = '';
             D.normalize = [];
                     D.p = cell(0);
                  D.path = pwd;
                 D.pcond = cell(0);
                 D.prand = cell(0);
                    D.pu = cell(0);
                    D.py = cell(0);
                D.py_sep = struct;
                 D.pystd = cell(0);
             D.rand_type = [];
    D.response_parameter = '';
                     D.t = '';
                  D.tExp = [];
                  D.tLim = [];
               D.tLimExp = [];
                D.tUnits = cell(0);
                     D.u = cell(0);
                D.uNames = cell(0);
                     D.y = cell(0);
                  D.yExp = [];
               D.yExpRaw = [];
               D.yExpStd = [];
            D.yExpStdRaw = [];
                D.yNames = cell(0);
                D.yUnits = cell(0);
%               D.checkstr = '';
%                    D.fkt = '';
%                  D.cLink = [];

                    
end