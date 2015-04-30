function str = arAssembleReactionStr(source, target, simbio)

if(~exist('simbio', 'var'))
    simbio = false;
end
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
if(simbio)
    if(isempty(s))
        s = 'null';
    end
    if(isempty(t))
        t = 'null';
    end
end
    
str = sprintf('%s -> %s', s, t);

        
