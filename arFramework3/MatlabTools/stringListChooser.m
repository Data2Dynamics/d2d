function out = stringListChooser(liste, default, zeigen)

if(~exist('default', 'var'))
    default = 1;
end
if(~exist('zeigen', 'var'))
    zeigen = true;
end

if(zeigen)
    maxlen = max(cellfun(@length,liste));
    anno = cell(1,length(liste));
    for j=1:length(liste)
        % the headers are the same and should be generated at the same place where anno is generated, overwriting is no problem
        [anno{j},vals{j},header] = readParameterAnnotation(liste{j});  % vals contains npara, nfitted, chi2, ... and is useful when working with a breakpoint
    end        
     
    maxlen2 = max(cellfun(@length,anno));
    fprintf(['         %',sprintf('%i',maxlen),'s%s\n'],' ', header(1:min(length(header),maxlen2)));
    for j=1:length(anno)
        fprintf(['#%3i : %-',num2str(maxlen),'s  %s\n'], j, liste{j}, anno{j});
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

function [anno,vals,header] = readParameterAnnotation(filename_tmp)
filename_pars = ['./Results/' filename_tmp '/workspace_pars_only.mat'];
vals = struct;

if(exist(filename_pars,'file'))
    S = load(filename_pars);
    
    if isfield(S.ar,'checkstrs')
        checkstr = S.ar.checkstrs;
    else
        checkstr = struct;
        checkstr.data = '';
        checkstr.para = '';
        checkstr.fitting = '';
        checkstr.fkt = '';
    end
    
    nstrH = '   N ';
    if(isfield(S.ar,'ndata'))
        nstr = sprintf('%4i  ',S.ar.ndata);
        vals.N = S.ar.ndata;
    else
        nstr = sprintf('%4s  ','NA');
        vals.N = NaN;
    end
    priorstrH = '#prior ';
    if(isfield(S.ar,'nprior'))
        priorstr = sprintf('%6i ',S.ar.nprior);
    else
        priorstr = sprintf('%6s ','NA');        
    end
    pstrH = '   #p ';
    if(isfield(S.ar,'p'))
        pstr = sprintf('%4i ',length(S.ar.p));
        vals.np = length(S.ar.p);
    else
        pstr = sprintf('%4s ','NA');
        vals.np = NaN;
    end
    qstrH = '#fitted ';
    if(isfield(S.ar,'qFit'))
        qstr = sprintf('%6i ',sum(S.ar.qFit==1));
        vals.nfit = sum(S.ar.qFit==1);
    else
        qstr = sprintf('%6s ','NA');
        vals.nfit = NaN;
    end
    
    errstrH = 'fitterrors';
    if isfield(S.ar,'config') && isfield(S.ar.config,'fiterrors')        
        errstr = sprintf('%9i',S.ar.config.fiterrors);
        if isfield(S.ar.config,'useFitErrorCorrection') && S.ar.config.useFitErrorCorrection
            errstr(1:6) = ' BessC';
        end
    else
        errstr = sprintf('%9s','NA');
    end

    chi2strH = sprintf('       Merit-fkt ');
    chi2str  = sprintf('       %9s  ','NA');
    if(isfield(S.ar,'chi2fit'))
        if(isfield(S.ar,'config') && S.ar.config.fiterrors == 1)
            chi2str = sprintf('-2*LL=%10g  ', ...
                2*S.ar.ndata*log(sqrt(2*pi)) + S.ar.chi2fit);
            vals.chi2 = 2*S.ar.ndata*log(sqrt(2*pi)) + S.ar.chi2fit;
        else
            chi2str = sprintf(' chi2=%10g  ', S.ar.chi2fit);
            vals.chi2 =  S.ar.chi2fit;
        end
    else
        vals.chi2 = NaN;
    end
    
    lhsstrH = ' #LHS';
    lhsstr = sprintf('%6s','NA');
    if(isfield(S.ar,'ps'))
        if ~isempty(S.ar.ps)
            lhsstr = sprintf('%6i',size(S.ar.ps,1));
        else
            lhsstr = sprintf('%6i',0);
        end
    end
    
    plestrH = ' #PLE';
    plestr = sprintf('%5i',0); % default
    if(isfield(S.ar,'ple'))
        if(isfield(S.ar.ple,'chi2s'))
            nple = sum(~cellfun(@isempty,S.ar.ple.chi2s));
            plestr = sprintf('%5i',nple);
        end
    elseif exist(['./Results/',filename_tmp,filesep,'workspace_pars_only.mat'],'file')==2
        % %% too slow:
        % vars = who('-file',['./Results/',filename_tmp,filesep,'workspace.mat'],'ple*');
        % if ~isempty(intersect(vars,'pleGlobals'))
        warning('off','MATLAB:load:variableNotFound'); % faster
        set(0, 'DefaultFigureVisible', 'off') % required for not displaying figures in the workspace, I found no nicer solution
        tmp = load(['./Results/',filename_tmp,filesep,'workspace_pars_only.mat'],'-mat','pleGlobals');
        set(0, 'DefaultFigureVisible', 'on') % required for not displaying figures in the workspace, I found no nicer solution
        warning('on','MATLAB:load:variableNotFound');
        if isfield(tmp,'pleGlobals');
            if ~isempty(tmp.pleGlobals)
                % disp('Old PLE as variable pleGlobals available. Transferred to ar.ple ...');
                nple = sum(~cellfun(@isempty,tmp.pleGlobals.chi2s));
                plestr = sprintf('(%3i)',nple);
            end
        end
    end
    
    anno = sprintf('%s%s%s%s%s%s%s%s  %s  %s  %s  %s',chi2str,nstr,pstr,qstr,priorstr,errstr,lhsstr,plestr,checkstr.para, checkstr.data, checkstr.fitting, checkstr.fkt);
    header = sprintf('%s%s%s%s%s%s%s%s%34s%34s%34s%34s',chi2strH,nstrH,pstrH,qstrH,priorstrH,errstrH,lhsstrH,plestrH,'parameter-settings', 'data-settings', 'fit-setting', 'setup-setting');
else
    anno = '';
    header = '';
end
