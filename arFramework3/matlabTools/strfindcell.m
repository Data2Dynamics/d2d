function q = strfindcell(data, pattern)

q = ~cellfun(@isempty, cellfun(@strfind, data, repmat({pattern},size(data,1),size(data,2)), 'UniformOutput', false));
