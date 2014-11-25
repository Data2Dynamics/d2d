function out = stringListChooser(liste, default, zeigen)

if(~exist('default', 'var'))
    default = 1;
end
if(~exist('zeigen', 'var'))
    zeigen = true;
end

if(zeigen)
    maxlen = max(cellfun(@length,liste));
    
    for j=1:length(liste)
        anno = readParameterAnnotation(liste{j});
        fprintf(['#%3i : %-',num2str(maxlen),'s  %s\n'], j, liste{j}, anno);
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

function anno = readParameterAnnotation(filename_tmp)
filename_pars = ['./Results/' filename_tmp '/workspace_pars_only.mat'];

if(exist(filename_pars,'file'))    
    S = load(filename_pars);

    nstr = '';
    priorstr = '';
    pstr = '';
    qstr = '';
    chi2str = '';
    errstr = '';
    
    if(isfield(S.ar,'ndata'))
        nstr = ['N=',sprintf('%4i ',S.ar.ndata),' '];
    end
    if(isfield(S.ar,'nprior'))
        priorstr = ['#prior=',sprintf('%3i ',S.ar.nprior),' '];
    end
    if(isfield(S.ar,'p'))
        pstr = ['#p=',sprintf('%3i ',length(S.ar.p)),' '];
    end
    if(isfield(S.ar,'qFit'))
        qstr = ['#fitted=',sprintf('%3i ',sum(S.ar.qFit==1)),' '];
    end
    if(isfield(S.ar,'chi2fit'))
        chi2str = ['chi2fit=',sprintf('%4.3f ',S.ar.chi2fit),' '];
    end    
    if(isfield(S.ar,'config'))
        if(S.ar.config.fiterrors==1)
            errstr = 'errors fitted';
        end
    end
    
    anno = sprintf('(%s%s%s%s%s%s)',chi2str,nstr,pstr,qstr,priorstr,errstr);
else
    anno = '';
end
 

 