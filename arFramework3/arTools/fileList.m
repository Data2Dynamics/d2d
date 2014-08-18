% Dateiliste erstellen

function out = fileList(filepath, searchpattern, doAnd)

if(~exist('searchpattern', 'var'))
    searchpattern = {};
end
if(~exist('doAnd', 'var'))
    doAnd = false;
end

filesyslist = dir(filepath);
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