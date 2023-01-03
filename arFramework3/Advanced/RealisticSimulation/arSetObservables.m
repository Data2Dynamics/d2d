function arSetObservables(cond,obslog)
%OBSERVABLES Define Observables
%   Detailed explanation goes here
if ~exist('cond','var') || isempty(cond)
    cond=1;
end
if ~exist('obslog','var') || isempty(obslog)
    fprintf('No information about logging given. Logarithm of observables is drawn realisticly. \nWatch out! Logarithm of intensities can cause differences to original model. If this is not what you want, set arRealisticDesign([],false).\n')
    obslog = false;
end

global ar

data = ar.model(1).condition(cond).xFineSimu;
label = ar.model(1).x;
load([fileparts(mfilename('fullpath')) filesep 'ObservableDraw.mat']);

% Are there constant states ? Throw them out.
Const = [];
for i=1:size(data,2)
    if (abs(range(data(:,i))/max(data(:,i))) < 1e-8 || isinf(range(data(:,i))/max(data(:,i))) || isnan(range(data(:,i))/max(data(:,i))) || abs(range(data(:,i)))<1e-8 )
        Const = [Const i];
    end
end
data(:,Const) = [];
label(Const) = [];
ns = length(label);
if ns == 0
    error('No dynamics found, just constants. Check in arSetObservables.m')
end
if ns<5
    error('Realistic Design not trained for less than 5 dynamic species.')
end

% Draw number of observables
    nobs = round(obs.obs(randi([1 length(obs.obs)]))*ns);
    if nobs>ns
        nobs = ns;
    end

    % compounds
    r = randi([1 length(obs.comp)]);
    nc = round(obs.comp(r)*nobs);
    ncadd = round(obs.compadd(r));

    if nc>nobs
        nc=nobs;
    end
    if ncadd>ns
        ncadd = ns;
    end
    if nc>0
        [~,idcom] = sort(rand(ncadd,ns),2);
        idcom = sort(idcom(:,1:ncadd),2);
        idcom = unique(idcom,'rows');
        if size(idcom,1)<nc
            [~,idcom2] = sort(rand(ncadd,ns),2);
            idcom2 = sort(idcom2(:,1:ncadd),2);
            idcom = unique([idcom; idcom2],'rows');
            if size(idcom,1)<nc
                nc = size(idcom,1);
            end
        end
        if size(idcom,1)>nc
            idcom = idcom(randperm(size(idcom,1),nc),:);
        end
        ycom = zeros(size(data,1),nc);
        for icadd=1:ncadd
            ycom = ycom + data(:,idcom(:,icadd));
        end
        labelcomind = label(idcom);
    end

    % Log?
    nlog = obs.logobs(randi([1 length(obs.obs)]));
    if ~obslog
        nlog = 0;
    end
    if nlog>0.5
        nsca = round(obs.scalog(randi([1 length(obs.scalog)]))*nobs);
        noff = round(obs.offlogrel(randi([1 length(obs.offlogrel)]))*nsca);
        ninit = round(obs.initlog(randi([1 length(obs.initlog)]))*nobs);
    else
        nsca = round(obs.scanonlog(randi([1 length(obs.scanonlog)]))*nobs);
        noff = round(obs.offnonlogrel(randi([1 length(obs.offnonlogrel)]))*nsca);
        ninit = round(obs.initnonlog(randi([1 length(obs.initnonlog)]))*nobs);
    end


% indeces
idobs = randperm(ns,nobs-nc);
labelobs = label(idobs);
idlog = binornd(1,nlog,nobs,1);
idsca = randperm(nobs,nsca);
idoff = idsca(randperm(nsca,noff));
idinit = randperm(nobs,ninit);

if nc>0
    newdata = [data(:,idobs) ycom];
else
    newdata = data(:,idobs);
end

% Write obs in RealisticData.def
if ~exist(['.' filesep 'Data'],'dir')
    mkdir('Data')
end
fileID = fopen(['Data' filesep 'RealisticData.def'],'w'); 
fprintf(fileID,'%s\n','DESCRIPTION') ;
fprintf(fileID,'\n%s\n','PREDICTOR') ;
fprintf(fileID,'\n%s\t%s\t%s\t%s\t%i\t%i\n','time','T','n/a','time',0,ar.model.tLim(2)) ;
fprintf(fileID,'\n%s\n','INPUTS') ;
fprintf(fileID,'\n%s\n','OBSERVABLES') ;
% compounds
ylabel = labelobs;
for i=1:nc
    labelcom = labelcomind{i,1};
    comform = labelcomind{i,1};
    for ii=2:size(labelcomind,2)
        if ~isempty(labelcomind{i,ii})
            labelcom = [labelcom '_add_' labelcomind{i,ii}];
            comform = [comform '+' labelcomind{i,ii}];
        end
    end
    if length(labelcom)>63
        labelcom = labelcom(1:63);
    end
    ylabel = [ylabel labelcom];
    if idlog(i)
        logflag=1;
    else
        logflag=0;
    end
    if any(idsca==i)
        if any(idoff==i)
            fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelcom '_obs'],'C','n/a','conc.',0,logflag,['"offset_'  labelcom '+ scale_' labelcom '*( ' comform ')"']) ;
        else
            fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelcom '_obs'],'C','n/a','conc.',0,logflag,['"scale_' labelcom '*(' comform ')"']) ;
        end
    else
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelcom '_obs'],'C','n/a','conc.',0,logflag,['"' comform '"']) ;
    end
end
for i=1:nobs-nc
    if idlog(i+nc)
        logflag=1;
    else
        logflag=0;
    end
    if any(idsca==i+nc)
        if any(idoff==i+nc)
            fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelobs{i} '_obs'],'C','n/a','conc.',0,logflag,['"offset_' labelobs{i} '+ scale_' labelobs{i} '* ' labelobs{i} '"']) ;
        else
            fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelobs{i} '_obs'],'C','n/a','conc.',0,logflag,['"scale_' labelobs{i} '*' labelobs{i} '"']) ;
        end
    else
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelobs{i} '_obs'],'C','n/a','conc.',0,logflag,['"' labelobs{i} '"']) ;
    end
end

fprintf(fileID,'\n%s\n','ERRORS') ;
for i = 1:length(ylabel)
    fprintf(fileID,'%s\t%s\n',[ylabel{i} '_obs'],[ '"sd_' ylabel{i} '"']) ;
end
% Rel + abs
% for i = 1:length(yNames)
%     fprintf(fileID,'%s\t%s\n',yNames{i},[ '"sd_abs_' yNames{i} ' + sd_rel_' yNames{i} ' * ' yNames{i} '"']) ;
% end

fprintf(fileID,'\n%s\n','CONDITIONS') ;
fprintf(fileID,'\n%s\n','RANDOM') ;
fprintf(fileID,'\n%s\n','PARAMETERS') ;

for i=1:length(idinit)
    value = log10(nanmin(newdata(newdata(:,idinit(i))>0,idinit(i)))/2);
    fprintf(fileID,'%s\t%f\t%i\t%i\t%f\t%f\n',['init_' ylabel{idinit(i)}],value,1,1,floor(value-2),ceil(value+2));
end

% Error parameters
sdmodel = -0.96 + randn(1)*0.3;
sd = nan(nobs,1);
for i = 1:nobs
    sd(i) = sdmodel + randn(1)*0.014;
    if ~idlog(i)
        sd(i) = sd(i) + log10(nanmean(newdata(:,i)));
    end
    fprintf(fileID,'%s\t%f\t%i\t%i\t%f\t%f\n',['sd_' ylabel{i}],sd(i),1,1,floor(sd(i)-2),ceil(sd(i)+2));
        %fprintf(fileID,'%s\t%f\t%i\t%i\t%f\t%f\n',['sd_' ylabel{i}],sdmean+randn(1)*0.4,1,1,-5,3);
end
fclose(fileID);

% Write in ar
ar.model.data.yFineSimu = newdata;
ar.model.data.y = ylabel; 
ar.model.data.tFine = ar.model.condition.tFine;

fprintf('Observables assigned. \n');



if naive
% kein scale/offset faktoren, keine compounds (just take first one)
% init_parameter same, error parameter same

    % Write obs in RealisticData_naive.def
    fileID = fopen(['Data' filesep 'RealisticData_naive.def'],'w'); 
    fprintf(fileID,'%s\n','DESCRIPTION') ;
    fprintf(fileID,'\n%s\n','PREDICTOR') ;
    fprintf(fileID,'\n%s\t%s\t%s\t%s\t%i\t%i\n','time','T','n/a','time',0,ar.model.tLim(2)) ;
    fprintf(fileID,'\n%s\n','INPUTS') ;
    fprintf(fileID,'\n%s\n','OBSERVABLES') ;
    % compounds
    ylabel = labelobs;
    for i=1:nc
        labelcom = labelcomind{i,1};
        comform = labelcomind{i,1};
        ylabel = [ylabel labelcom];
        if idlog(i)
            logflag=1;
        else
            logflag=0;
        end
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelcom '_obs'],'C','n/a','conc.',0,logflag,['"' comform '"']) ;
    end
    for i=1:nobs-nc
        if idlog(i+nc)
            logflag=1;
        else
            logflag=0;
        end
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', [labelobs{i} '_obs'],'C','n/a','conc.',0,logflag,['"' labelobs{i} '"']) ;
    end

    fprintf(fileID,'\n%s\n','ERRORS') ;
    for i = 1:length(ylabel)
        fprintf(fileID,'%s\t%s\n',[ylabel{i} '_obs'],[ '"sd_' ylabel{i} '"']) ;
    end
    % Rel + abs
    % for i = 1:length(yNames)
    %     fprintf(fileID,'%s\t%s\n',yNames{i},[ '"sd_abs_' yNames{i} ' + sd_rel_' yNames{i} ' * ' yNames{i} '"']) ;
    % end

    fprintf(fileID,'\n%s\n','CONDITIONS') ;
    fprintf(fileID,'\n%s\n','RANDOM') ;
    fprintf(fileID,'\n%s\n','PARAMETERS') ;

    for i=1:length(idinit)
        value = log10(nanmin(newdata(newdata(:,idinit(i))>0,idinit(i)))/2);
        fprintf(fileID,'%s\t%f\t%i\t%i\t%f\t%f\n',['init_' ylabel{idinit(i)}],value,1,1,floor(value-2),ceil(value+2));
    end

    % Error parameters
%     for i = 1:nobs
%        fprintf(fileID,'%s\t%f\t%i\t%i\t%f\t%f\n',['sd_' ylabel{i}],sd(i),1,1,floor(sd(i)-2),ceil(sd(i)+2));
%     end
    fclose(fileID);
    
    fprintf('Observables without scaling or compound saved to "RealisticData_naive.def". \n');
end