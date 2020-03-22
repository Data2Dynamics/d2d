% status = arLoadLatest([pattern])
% 
%  Loads the most recent workspace. With the regular expression "pattern"
%  the most recent workspace which fits to the regular experession
%  specified by pattern is loaded.
% 
%  The date of the created folder ar.config.savepath (i.e. in the Results
%  folder) is evaluated (if the workspace.mat is overwritten later, this
%  does not count in the current implementation).
% 
%     pattern          optional regular experession 
% 
%     status           optional boolean, specifies if a workspace was found and loaded 
%
% Examples:
% 
%   arLoadLatest
% 
%   arLoadLatest('PLE');

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
datNum = [d.datenum];
dates = dates([d.isdir]);
datNum = datNum([d.isdir]);
folders = names([d.isdir]);

folders = folders(3:end); % remove . and ..
dates = dates(3:end); % remove . and ..
datNum = datNum(3:end); % remove . and ..

% only folders which contain workspace.mat
ok = false(size(folders));
for i=1:length(folders)
    ok(i) = exist(['Results',filesep,folders{i},filesep,'workspace.mat'],'file');
end

folders = folders(ok);
dates = dates(ok);
datNum = datNum(ok);

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
    datNum = datNum(ind);
end

if ~isempty(dates)
    [~,rf] = sort(datNum); % sorting according to date  
    arLoad(folders{rf(end)});
    status = true;
else
    status = false;
    disp('No proper workspace found.')
end

if nargout>0
    varargout{1} = status;
end
