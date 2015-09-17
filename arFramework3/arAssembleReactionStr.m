function str = arAssembleReactionStr(source, target, simbio, sourceCoeffs, targetCoeffs)

if(~exist('simbio', 'var'))
    simbio = false;
end
if(isempty(source))
    s = '';
else
    s = source{1};
end
if ( nargin > 3 )
    for jj=2:length(source)
        s = [s ' + ' num2str(sourceCoeffs(jj)) + ' ' + source{jj}]; %#ok<AGROW>
    end
else
    for jj=2:length(source)
        s = [s ' + ' source{jj}]; %#ok<AGROW>
    end    
end
if(isempty(target))
    t = '';
else
    t = target{1};
end
if ( nargin > 4 )
    for jj=2:length(target)
        t = [t ' + ' num2str(targetCoeffs(jj)) + ' ' + target{jj}]; %#ok<AGROW>
    end
else
    for jj=2:length(target)
        t = [t ' + ' target{jj}]; %#ok<AGROW>
    end    
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

        
