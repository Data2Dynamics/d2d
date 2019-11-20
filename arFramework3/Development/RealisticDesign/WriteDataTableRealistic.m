function [T, yExp] = WriteDataTableRealistic( tT, y, yNames, folder_name )
%WRITEDATATABLE Summary of this function goes here
%   Detailed explanation goes here

global ar

if(~exist('tT', 'var') || isempty(tT))
    error('No TimePoints in WriteDataTable function')
end
if(~exist('yNames', 'var') || isempty(yNames))
    error('No Observable Names in WriteDataTable function')
end
if(~exist('y', 'var') || isempty(y))
    error('No simulated ObsData in WriteDataTable function')
end
% if(~exist('dt', 'var') || isempty(dt))
%     dt = zeros(length(yNames),1);
% %     error('No time spacing of deleted first steady points in WriteDataTable function')
% end
if(~exist('folder_name', 'var') || isempty(folder_name))
    folder_name = '';
end


% Replace too small values
y(y<1e-8 & y>0)=1e-8;
y(y>-1e-6 & y<0)=1e-8;
if length(y(y<1e-7))<size(y,2)
    y(y<1e-7)=1e-7;
end
if length(y(y<1e-6))<size(y,2)
    y(y<1e-6)=1e-6;
end
if length(y(y<1e-5))<size(y,2)
    y(y<1e-5)=1e-5;
end
if length(y(y<1e-4))<size(y,2)
    y(y<1e-4)=1e-4;
end
if length(y(y<1e-3))<size(y,2)
    y(y<1e-3)=1e-3;
end
% for i=1:size(y,2)
%     if any(y(:,i)<0) && abs(min(y(:,i)))<max(y(:,i))/2
%         y(:,i) = abs(y(:,i));
%     end
% end

% Order data all in one matrix
T = unique(sort(tT(~isnan(tT(:)))));
yExp = nan(size(T,1),size(tT,2));
for i=1:size(y,1)
    for j=1:size(y,2)
        for k=1:size(T,1)
            if tT(i,j)==T(k)
                yExp(k,j) = y(i,j);
            end
        end
    end
end
while all(isnan(yExp(end,:)))
    yExp(end,:) = [];
    T(end) = [];
end

Data = nan(size(T,1),size(yExp,2)+1);
Text = cell(1,length(yNames)+1);

Text{1,1} = 't';
Data(:,1) = T;

Text(1,2:size(yExp,2)+1) = yNames(:,1);
Data(:,2:size(yExp,2)+1) = yExp;


Raw = [Text; num2cell(Data)];
xlswrite(['RealisticDesign' folder_name '/RealisticData.xls'],Raw);

% Save in struct
ar.model.data.tExp = T;
ar.model.data.yExp = yExp;
% Set first data point at init
%ar.model.data.yExp(1,:) = ar.yinit;


end

