%   arAddDataStruct(D, [m])
% 
% This function adds a data struct and appends it to ar.model(m)
% This function might serve as a replacement or arLoadData. Compiling has
% to be performed after adding the/all data structs.
% 
% 
%   D       data struct like usually found in ar.model.data
%   m       model index
%           Default: m = length(ar.model) [like in arLoadData]
% 

function arAddDataStruct(D, m)
global ar

if ~exist('m','var') || isempty(m)
    m = length(ar.model);
end

if isfield(ar.model(m),'data')
    ar.model(m).data(end+1) = D;
else 
    ar.model(m).data = D;
end

%% Create the respective plot-struct:
d = length(ar.model(m).data);
if(~isfield(ar.model(m), 'plot'))
    ar.model(m).plot(1).name = D.name;
else
    ar.model(m).plot(end+1).name = D.name;
end
ar.model(m).plot(end).doseresponse = ar.model(m).data(d).doseresponse;
ar.model(m).plot(end).doseresponselog10xaxis = true;
ar.model(m).plot(end).dLink = d;
ar.model(m).plot(end).ny = length(ar.model(m).data(d).y);
ar.model(m).plot(end).condition = {};


% remember the function call
ar.setup.commands{end+1} = mfilename; % this file name
ar.setup.arguments{end+1} = {D,m}; % 
ar.setup.datafiles{end+1} = '';
ar.setup.modelfiles{end+1} = '';


