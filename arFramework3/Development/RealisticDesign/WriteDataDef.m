
function WriteDataDef(plog,p2)

global ar

yNames = ar.model.data.y;
tmax = max(nanmax(ar.model.data.tExp));

if ~exist('tmax','var') || isempty('tmax')
    error('tmax has to be the same as in Biomodel, to ensure the right inputs in the simulation: yFine, tFine.')
end    
if ~exist('plog','var')
    plog = false;               % Log Parameters ? Just if you sure it works.
end    
if ~exist('p2','var')
    p2 = false;               % Log Parameters ? Just if you sure it works.
end   
fclose('all');
% delete pre-existing file
if exist([pwd 'Data/RealisticData.def'],'file')
    delete(['Data/RealisticData.def'])
end
if exist([pwd 'Data/RealisticData.xls'],'file')
    delete(['Data/RealisticData.xls'])
end

% Write Data .def with new observables
fileID = fopen(['Data/RealisticData.def'],'w'); 
fprintf(fileID,'%s\n','DESCRIPTION') ;
fprintf(fileID,'\n%s\n','PREDICTOR') ;
fprintf(fileID,'\n%s\t%s\t%s\t%s\t%i\t%i\n','time','T','n/a','time',0,tmax) ;
fprintf(fileID,'\n%s\n','INPUTS') ;
fprintf(fileID,'\n%s\n','OBSERVABLES') ;
c=1; scaled = cell(1);
for i = 1:length(yNames)
    % if scaled
    expression = 'scaled';
    replace = '';
    newStr = regexprep(yNames{i,1},expression,replace);
    if strcmp(yNames{i,1},newStr)==0
        expression_new = '_obs';
        replace_new = '';
        newStr_new = regexprep(newStr,expression_new,replace_new);
        expression_sl = ' ';
        replace_sl = '';
        newStr_sl = regexprep(newStr_new,expression_sl,replace_sl);
        scaled{c} = newStr_sl;
        c=c+1;
        yNames{i} = [newStr_sl '_obs'];
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n',yNames{i},'C','n/a','conc.',0,0,['"offset_' newStr_sl ' + scale_' newStr_sl '*' newStr_sl '"']) ;
    end
    % if additive compound
    expression_a = '+';
    replace_a = '_';
    newStr_a = strsplit(yNames{i,1},'+');
    if length(newStr_a)>1
        expression_new = '_obs';
        replace_new = '';
        newStr_new = regexprep(newStr_a,expression_new,replace_new);
        expression_al = ' ';
        replace_al = '';
        newStr_al = regexprep(newStr_new,expression_al,replace_al);
        yNames{i} = [newStr_al{1} '_add_' newStr_al{2} '_obs'];
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n',yNames{i},'C','n/a','conc.',0,0,[ '"' newStr_al{1} '+' newStr_al{2} '"']) ;
    end
    
    % Write rest yNames in def
    if strcmp(yNames{i,1},newStr)~=0 && length(newStr_a)<=1
        expression = '_obs';
        replace = '';
        newStr = regexprep(yNames{i,1},expression,replace);
        fprintf(fileID,'%s\t%s\t%s\t%s\t%i\t%i\t%s\n', yNames{i,1},'C','n/a','conc.',0,0,[ '"' newStr '"']) ;
    end
end

fprintf(fileID,'\n%s\n','ERRORS') ;
% Abs
for i = 1:length(yNames)
    fprintf(fileID,'%s\t%s\n',yNames{i,1},[ '"sd_rel_' yNames{i,1} ' * ' yNames{i,1} ' + 1e-7"']) ;
end
% Rel + abs
% for i = 1:length(yNames)
%     fprintf(fileID,'%s\t%s\n',yNames{i,1},[ '"sd_abs_' yNames{i,1} ' + sd_rel_' yNames{i,1} ' * ' yNames{i,1} '"']) ;
% end

fprintf(fileID,'\n%s\n','CONDITIONS') ;
if ~isempty(scaled)
    for i=1:length(scaled) 
        fprintf(fileID,'%s\t%s\n',['init_' scaled{i}], ['"init_' scaled{i} ' * scale_' scaled{i} '"']);
    end
end

fclose(fileID);


%% Write model .def
% for parameter changes (set bounds, log them)

file = dir('Models/BIOMD*');
fin = fopen(['Models/' file(1).name]);
s = textscan(fin, '%s','delimiter','\n');
idx = find(strcmp(s{1}, 'PARAMETERS'), 1, 'last');
if isempty(idx)
    idx=0;
end
clear s
fclose(fin);

name = 'Realistic';
if plog
    name = [name '_Log'];
end
if p2
    name = [name '_p2'];
end
fout = fopen(['Models/' name '.def'], 'wt');
fin = fopen(['Models/' file(1).name]);
tline = fgetl(fin);

count = 0;
while ischar(tline) && count < idx 
    if ischar(tline)
        fprintf(fout, '%s\n', tline);
    end
    tline = fgetl(fin);
    count = count + 1;
end

% Get parameters from Biomodel, (log them), set bounds p-2<p<p+2
% PARAMETERS
% pExternLabels   pExtern    qFitExtern    qLog10Extern    lbExtern    ubExtern
% MinMagn = round(log10(min(min(ar.model.data.yFineSimu(ar.model.data.yFineSimu>0))))-1,1,'significant');
MinMagn = round(min(ar.Magn),1,'significant')-2;
if p2
    if plog
        for i=1:length(ar.pExtern)
                if abs(ar.pExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},MinMagn,1,1,MinMagn-2,MinMagn+2);
                else
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},log10(ar.pExtern(i)),1,1,log10(ar.pExtern(i))-2,log10(ar.pExtern(i))+2);
                end
        end
    else
        for i=1:length(ar.pExtern)
                if abs(ar.pExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},10.^(MinMagn),1,0,10.^(MinMagn),10.^(MinMagn+4));
                else
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},ar.pExtern(i),1,0,ar.pExtern(i)/100,ar.pExtern(i)*100);
                end
        end
    end
else
    if plog
        for i=1:length(ar.pExtern)
                if abs(ar.pExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},MinMagn,1,1,MinMagn,log10(ar.ubExtern(i)));
                elseif abs(ar.lbExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},log10(ar.pExtern(i)),1,1,MinMagn,log10(ar.ubExtern(i)));     
                else
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},log10(ar.pExtern(i)),1,1,log10(ar.lbExtern(i)),log10(ar.ubExtern(i)));
                end
        end
    else
        for i=1:length(ar.pExtern)
                if abs(ar.pExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},10.^(MinMagn),1,0,10.^(MinMagn),ar.ubExtern(i));
                elseif abs(ar.lbExtern(i))<10^(MinMagn)
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},ar.pExtern(i),1,0,min(10.^(MinMagn),ar.pExtern(i)),ar.ubExtern(i));
                else
                    fprintf(fout,'\n%s\t%.3g\t%i\t%i\t%.3g\t%.3g',ar.pExternLabels{i},ar.pExtern(i),1,0,ar.lbExtern(i),ar.ubExtern(i));
                end
        end
    end
end
fclose(fin);
fclose(fout);




% %% Compile
% % has to be! Because error model of sd * Obs has to be implemented before
% % arSimuData
% 
% if exist('Results/Realistic','dir')
%     rmdir('Results/Realistic','s')
% end
% fprintf('Compile new Observables and Parameters.')
% arInit
% ar.config.checkForNegFluxes = false;
% name = 'Realistic';
% if plog
%     name = [name '_Log'];
% end
% if p2
%     name = [name '_p2'];
% end
% if exist(['Models/Realistic' name '.def'],'file')
%     arLoadModel('Realistic')
% else
%     error('Couldn''t find model.def in WriteDataDef.m');
% end
% arLoadData(['RealisticData']);
% arCompileAll(true);
% 
% % save
% for i=1:length(ar.p)
%     if strncmp(ar.pLabel{i},'offset_',6)
%         ar.lb(i) = ar.p(i);
%         ar.ub(i) = ar.p(i)+2;
%     elseif strncmp(ar.pLabel{i},'scale_',6)
%         ar.qLog10(i) = 1;
%         ar.lb(i) = -1;
%         ar.p(i) = 0;
%         ar.ub(i) = 1;
%     elseif strncmp(ar.pLabel{i},'sd_',3)
%         ar.lb(i) = -1.2;
%         ar.ub(i) = -0.8;
%     end
% end
%     
% arSave(['Realistic' name]);
% [~,ws]=fileparts(ar.config.savepath);
% movefile(['Results/' ws],['Results\Realistic' name]);
% fprintf('Offset and scaling factors are set.\n');
% fprintf(['Realistic workspace for model ' file(1).name(1:end-4) ' saved to file ./Results/Realistic' name '/workspace.mat \n']);
