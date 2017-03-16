% Dateiliste erstellen

function out = fileList(filepath, searchpattern, doAnd, onlyFolders)

if(~exist('searchpattern', 'var')  || isempty(searchpattern))
    searchpattern = {};
end
if(~exist('doAnd', 'var') || isempty(doAnd))
    doAnd = false;
end
if(~exist('onlyFolders', 'var') || isempty(onlyFolders))
    onlyFolders = false;
end

if(~iscell(searchpattern))
    searchpattern = {searchpattern};
end

filesyslist = dir(filepath);
if onlyFolders
    filesyslist = filesyslist([filesyslist.isdir]==1);  % only show folders
end
out = {};
count = 0;

for j=1:length(filesyslist)
    q = true;
    if(length(searchpattern)>0) %#ok<ISMT>
        q = findPattern(filesyslist(j).name, searchpattern, doAnd);
    end
    if(filesyslist(j).name(1)~='.' && q)
        count = count + 1;
        out{count} = filesyslist(j).name; %#ok<AGROW>
    end
end

%% place "BestFit" always on Top (in order to have sorting according to date)
[~,~,ia] = intersect(lower('BestFit'),lower(out));
if ~isempty(ia)
    out = [out(ia),out(setdiff(1:length(out),ia))];
end



function out = findPattern(string, searchpattern, doAnd)

count = 0;
out = false;

for j=1:length(searchpattern)
    count = count + length(strfind(string, searchpattern{j}));
end

if(~doAnd)
    if(count>0)
        out = true;
    end
else
    if(count==length(searchpattern))
        out = true;
    end
end