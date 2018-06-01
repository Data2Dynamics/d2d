% arLoadLatest
% 
%  Loads the most recent workspace 
% 
% arLoadLatest(pattern)
% 
%   loads the most recent workspace which fits to the regular experession
%   specified by pattern
% 
% 
% 
% Examples:
% 
% arLoadLatest
% 
% arLoadLatest('PLE');

function arLoadLatest(pattern)
if ~exist('pattern','var') || isempty(pattern)
    pattern = [];
end

d = dir('Results');
names = {d.name};
dates = {d.date};
dates = dates([d.isdir]);
folders = names([d.isdir]);

if ~isempty(pattern)
    ind = find(~cellfun(@isempty,regexp(folders,pattern)));
    if isempty(ind)
        warning('Pattern ''%s'' does not match to any workspace in folder Results => nothing loaded.',pattern);
        return
    end
    folders = folders(ind);
    dates = dates(ind);
end

[~,rf] = sort(datenum(dates));

arLoad(folders{rf(end)});


