function str = arAssambleReactionStr(source, target)

if(isempty(source))
    s = '';
else
    s = source{1};
end
for jj=2:length(source)
    s = [s ' + ' source{jj}]; %#ok<AGROW>
end
if(isempty(target))
    t = '';
else
    t = target{1};
end
for jj=2:length(target)
    t = [t ' + ' target{jj}]; %#ok<AGROW>
end
str = sprintf('%s -> %s', s, t);
