%% SortData
% This file contains code to structure each data file and split it into one excel file for each dose. Additionally WT and KO are splitted into different columns

cd Data; % Change to data path
file = 'ELISA-Nigericin-formatted'; % file

%%  Load data, get headers and make subfolders
folder = ['Doses/',file,'-Doses/'];
mkdir(folder);
data = xlsread([file,'.xlsx']);

%% Load data, get headers and make subfolders
[d,s] = xlsread(file);
header = s(1,:);
index_dose = find(contains(header,'dose'));
index_Response = find(contains(header,'Response'));
index_isKO = find(contains(header,'isKO'));
index_Response2 = length(header)+1;

doses_len = length(unique(data(:,index_dose)));
doses_num = unique(data(:,index_dose));
doses_str = num2str(doses_num);

%% Make empty files for each dose with headers
for i = 1:doses_len
    % name files fo reach dose
    i_str = strtrim(doses_str(i,:));
    file_xls{i} = [folder,file,'-Dosis-', i_str, '.xlsx'];
    file_def{i} = [folder,file,'-Dosis-', i_str, '.def'];
    % plain .def files for each dose:
    copyfile([file,'.def'], file_def{i}); 
    
    T = readtable([file,'.xlsx']);
    T(1:end,:) = [];
    ind = find(strcmp(T.Properties.VariableNames, {'Response'}));
    T.Properties.VariableNames{ind} = 'Response1';
    T.Response2 = zeros(0); % Add empty column 'Response2' for KO
    writetable(T,file_xls{i}); 
end

%% Fill the excel files for each dose
for i=1:doses_len
    i_str = strtrim(doses_str(i,:));
    i_num = doses_num(i);
    % split by doses
    ind = data(:, index_dose) == i_num;
    data_dose = data(ind,:);
    % split by isKO in different columns
    ind = data_dose(:,index_isKO) == 1;
    data_dose(ind,index_Response2) = data_dose(ind,index_Response);
    data_dose(ind,index_Response) = NaN;
    data_dose(~ind,index_Response2) = NaN;
    % write in own xls for each dose
    writematrix(data_dose,file_xls{i},'Sheet',1,'Range','A2')
end
