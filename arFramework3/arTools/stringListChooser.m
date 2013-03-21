function out = stringListChooser(liste, default, zeigen)

if(~exist('default', 'var'))
    default = 1;
end
if(~exist('zeigen', 'var'))
    zeigen = true;
end

if(zeigen)
    for j=1:length(liste)
        fprintf('#%3i : %s\n', j, liste{j});
    end
end

out = -1;
if(default~=0) % force input
    while(out<1 || out>length(liste))
        out = numberChooser2(sprintf('Please choose (1-%i) ', length(liste)), default, liste);
        if(out<1 || out>length(liste))
            fprintf('The number has to be in the range of %i-%i!\n', 1, length(liste));
        end
    end
else % allow zero input
    while(out<0 || out>length(liste))
        out = numberChooser2(sprintf('Please choose (1-%i) ', length(liste)), default, liste);
        if(out<0 || out>length(liste))
            fprintf('The number has to be in the range of %i-%i!\n', 0, length(liste));
        end
    end
end


function out = numberChooser2(text, default, liste)

if(~exist('text', 'var'))
    text = '';
end

out = nan;

while(isnan(out))
    defaultstr = '';
    if(exist('default', 'var'))
        defaultstr = [text '[ENTER = ' sprintf('%i', default) ',  ls = list]'];
    end
    eingabe = input([defaultstr ': '], 's');
    if(isempty(eingabe))
        if(exist('default', 'var'))
            out = default;
        else
            fprintf('This is not a number!\n');
            out = nan;
        end
    else
        if(strcmp(eingabe, 'ls'))
            for j=1:length(liste)
                fprintf('#%3i : %s\n', j, liste{j});
            end
            out = nan;
        else
            out = str2double(eingabe);
            if(isempty(out) || isnan(out))
                fprintf('This is not a number!\n');
                out = nan;
            end
        end
    end
end
