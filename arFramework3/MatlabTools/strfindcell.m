function q = strfindcell(data, pattern)

if(iscell(pattern))
    q = false(size(data));
    for j=1:length(pattern)
        q = q | strfindcell(data, pattern{j});
    end
    return
end

q = ~cellfun(@isempty, cellfun(@strfind, data, repmat({pattern},size(data,1),size(data,2)), 'UniformOutput', false));
