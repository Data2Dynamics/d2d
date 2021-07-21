% D = arReadDataFile(name, [extension], [opts])
% 
% read data file and pass it out as an struct
% since there is no information about the model details in this function, data file is 
% just stored as raw data
%
%   name             char    data definition file name
%                            if it is empty, loading data ignored
%   extension        char    data file extension: 'xls', 'xlsx' or 'csv'
%                            'none' = don't load data
%                            ['xls']
%   opts             strucr  to extract data path
%                            ['Data/']

function D = arReadDataFile(name, extension, opts)

if ~exist('opts','var') || isempty(opts.datapath_args)
    DataPath = 'Data/';
else
    DataPath = opts.datapath_args;
    if DataPath(end)~='/' && DataPath(end)~='\'
        DataPath = [DataPath,'/'];
    end
end

if ~exist(DataPath,'dir')
    error('folder %s does not exist',DataPath)
end

if ~ischar(name)
    error('class(name)~=''char''');
end
if strcmp(strrep(name,' ',''),name)~=1
    disp(name)
    error('File names should not contain empty spaces. Please remove it.');
end

if ~ischar(name) || ~ischar(extension) || nargin>3
    error(['arReadDataDef(name, extension , opts) input argument is deprecated !!! ' ...
        'Please see new usage of arReadDataDef(name, extension, opts) and function help text.']);
end


% not loading data
if ~exist('name','var') || (exist('name','var') && isempty(name))
    if exist('extension','var') && strcmp(extension,'none')
        D.name = 'NoName';
        D.extension  = extension ;
        warning('No data file loaded because no data file name introduced (there is no file name specified) and extension set to ''none'' ');
    elseif exist('extension','var') && (strcmp(extension,'xls') || strcmp(extension,'xlsx') || strcmp(extension,'csv'))
        D.name = 'NoName';
        D.extension  = extension ;
        warning('No data file loaded because no data file name introduced (there is no file name specified), even though the extension set to ''%s'' ', extension);
    else
        D.name = 'NoName';
        D.extension  = 'none';
        warning('Ignore extension %s, data not loaded because no data file name introduced (there is no file name specified)', extension);
    end    
    D.path = [pwd,filesep,DataPath];
    D.xlstimevar = {};
    D.xlstimes = [];
    D.xlsheader = {};
    D.xlsdata = [];
    D.xlsdataCell = {};
    orderfields(D);
    return
else
    if exist('extension','var') && strcmp(extension,'none')
        D.name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
        D.extension  = extension ;        
        D.path = [pwd,filesep,DataPath];
        D.xlstimevar = {};
        D.xlstimes = [];
        D.xlsheader = {};
        D.xlsdata = [];
        D.xlsdataCell = {};
        orderfields(D);
        warning('Ignore loading data file');
        return        
    end  
end


if ~exist('extension','var') || isempty(extension)
    extension = 'xls';
    % auto-select extension if not specified
    if exist([DataPath, name '.xlsx'],'file')
        extension = 'xlsx';
    elseif exist([DataPath, name '.xls'],'file')
        extension = 'xls';
    elseif exist([DataPath, name '.csv'],'file')
        extension = 'csv';
    end
end

if ~exist([DataPath, name, '.', extension],'file')
    %arParsingError(fid, 'data file corresponding to %s does not exist in folder %s/', name, DataPath)
    error('Data file corresponding to %s.%s does not exist in folder %s', name, extension, DataPath)
end

fid = fopen([DataPath, name, '.' , extension], 'r');
arFprintf(3, 'Open Data file %s%s.%s [OK]\n', DataPath, name, extension);

D.name = strrep(strrep(strrep(strrep(name,'=','_'),'.',''),'-','_'),'/','_');
D.path = [pwd,filesep,DataPath];
D.extension = extension;



% READ XLS FILE -----------------------------------------------------------

arFprintf(2, 'Loading data from %s%s.%s...', DataPath, name, extension);

% read from file
if contains(extension,'xls')

    warning('off','all')
    
    arFprintf(3, '[ OK ]\nBegin reading data (xls) ...');  
    T = readtable([DataPath, name, '.', extension],'ReadRowNames',false);   
    arFprintf(3, '[ OK ]\n');
    
    header = T.Properties.VariableNames;
    data = T.Variables;
    
    timevar = header(1);
    header = header(2:end);
    header = strrep(header,' ','');
    times = data(:,1);
    data = data(:,2:end);
    
    % remove data without headers
    idxVar=[];
    for i=1:length(header)
        if startsWith(header{i},'Var')
            idxVar(end+1)=i;
        end
    end
    header(idxVar)=[];
    data(:,idxVar)=[];
       
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
        
elseif strcmp(extension,'csv')
    arFprintf(3, '[ OK ]\nBegin reading data (csv) ...');
    [header, data, dataCell] = arReadCSVHeaderFile([DataPath, name '.csv'], ',', true);
    arFprintf(3, '[ OK ]\n');
    
    timevar = strtrim(header(1));
    header = header(2:end);
    times = data(:,1);
    data = data(:,2:end);
    dataCell = dataCell(:,2:end);
end

if ismember(timevar, header)
    arParsingError( fid, 'predictor variable repetition in data file header')
else
    D.xlstimevar = timevar;
end

D.xlstimes = times;
for i=1:size(dataCell,2)
    q = ismember(header, header(i));
    if sum(q)==1
        D.xlsheader(i) = header(i);
        D.xlsdata(:,i) = data(:,i);
        D.xlsdataCell(:,i) = dataCell(:,i);
    else
        arParsingError(fid, 'multiple data colums for observable %s', header(i))
    end
end


arCheckReservedWords(D.xlsheader, 'data file header')

D = orderfields(D);

if ~isstruct(fid)
    fclose(fid);
end

arFprintf(2, '[OK]\n');

end


