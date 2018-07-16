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

function varargout = arLoadLatest(pattern)
if ~exist('pattern','var') || isempty(pattern)
    pattern = [];
end
if ~exist('Results','dir')
    error('No results folder exist. arLoadLatest can only be executed in a D2D working directory.')
end

d = dir('Results');
names = {d.name};
dates = {d.date};
dates = dates([d.isdir]);
folders = names([d.isdir]);

folders = folders(3:end); % remove . and ..
dates = dates(3:end); % remove . and ..

% only folders which contain workspace.mat
ok = false(size(folders));
for i=1:length(folders)
    ok(i) = exist(['Results',filesep,folders{i},filesep,'workspace.mat'],'file');
end

folders = folders(ok);
dates = dates(ok);

if ~isempty(pattern)
    ind = find(~cellfun(@isempty,regexp(folders,pattern)));
    if isempty(ind)
        warning('Pattern ''%s'' does not match to any workspace in folder Results => nothing loaded.',pattern);
        if nargout>0
            varargout{1} = false;
        end
        return
    end
    folders = folders(ind);
    dates = dates(ind);
end

if ~isempty(dates)
    [~,rf] = sort(datenum(dates)); % sorting according to date  
    arLoad(folders{rf(end)});
    status = true;
else
    status = false;
    disp('No proper workspace found.')
end

if nargout>0
    varargout{1} = status;
end
