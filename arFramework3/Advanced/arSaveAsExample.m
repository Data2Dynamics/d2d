%   arSaveAsExample
% 
%   arSaveAsExample(makePlots)
%         makePlots     Default: true (new plots are made for x, y, and profiles
%                       false: No (new) plots are generated and saved
% 
%  Function for copying the current model to arFramwork3/Examples
%  If your current model exhibts interesting features or problems, use this
%  function to copy it to the examples folder.
% 
%  The folder is added to .gitignore
% 
%  The following is copied:
%       ar.fkt
%       *.def and *.xls of the loaded models and data
%       d2d workspaces saved with arSave
%       plots

function arSaveAsExample(makePlots)
if ~exist('makePlots','var') || isempty(makePlots)
    makePlots = true;
end

global ar 

name = input('Model will be saved as Example. Please enter the name: ','s');
name = strrep(name,' ','_');
% arSave(name);
arSave('current');


if makePlots
    plePlotMulti([],true);
end
if isstruct(ar.ple) 
    if isfield(ar.ple,'ps') && ~isempty(ar.ple.ps)
        copyfile([ar.ple.savePath],[ar.config.savepath,filesep,'PLE'])
    end
end

if makePlots
    arQplot('xy');
    arPlot(true,true);
end


pfad = strrep(which('arSaveAsExample','-all'),'arSaveAsExample.m',['Examples',filesep,'_myOwnExamples']);
pfad0 = strrep(which('arSaveAsExample','-all'),'arSaveAsExample.m',['']);
if iscell(pfad) & length(pfad)>1
    fprintf('Several possible target folders:\n');
    for i=1:length(pfad)
        fprintf(' %i: %s\n',i,pfad{i});
    end
    
    ok = 0;
    while ok==0
        in = input('Please specify target (enter a number): ','s');
        try
            in = str2num(in);
            if isnumeric(in)
                ok = 1;
            else
                ok = 0;
            end
        catch
            ok = 0;
        end
    end
   
    pfad = pfad{in};
    pfad0 = pfad0{in};
elseif iscell(pfad)
    pfad = pfad{1};
    pfad0 = pfad0{1};
end

pfad = [pfad,filesep,name];

if ~exist([pfad0,filesep,'Examples',filesep,'_myOwnExamples'],'dir')
    mkdir([pfad0,filesep,'Examples',filesep,'_myOwnExamples'])
end

mkdir([pfad,filesep,'Models']);
mkdir([pfad,filesep,'Data']);
mkdir([pfad,filesep,'Results']);

copyfile(ar.config.savepath,[pfad,filesep,ar.config.savepath]);
copyfile([ar.fkt,'*'],[pfad]);


%% which data and model files:
ms = {ar.model.name};
ds = cell(size(ms));
for m=1:length(ar.model)    
    if isfield(ar.model(m),'data')
        if(isfield(ar.model(m).data,'name'))
            ds{m} = {ar.model(m).data.name};
            try % new code (unfortunately not tested, sorry)
                prands = [ar.model(m).data.prand];
                for ii=1:length(prands)
                    ds{m} = regexprep(ds{m},['_',prands{ii},'(\d)+'],'');
                end
            catch            % old code
                disp('Please check replacement code of random parameters in arSaveAsExample.m')
                [uni,ia,ib]= unique(regexprep(ds{m},'_nExpID(\d)+',''));
                ds{m} = uni(ib);  % replace zurueck
                ds{m} = ds{m}(sort(ia)); % nur die unique, aber in alter reihenfolge
            end
            
        else
            ds{m} = [];
        end
    end
end

for m=1:length(ms)
    copyfile(['Models',filesep,ms{m},'.def'], [pfad,filesep,'Models']);
    for d=1:length(ds{m})
        copyfile(['Data',filesep,ds{m}{d},'.*'], [pfad,filesep,'Data']);
    end
end


matlab_version = ver;
save([pfad,filesep,'matlab_version.mat'], 'matlab_version')

movefile([prefdir,filesep,'History.xml'],pfad);

%%
d = dir;
files = {d.name};
ind = find(~cellfun(@isempty,regexpi(files,'setup*.m')));
for i=1:length(ind)
    copyfile(files{ind(i)},pfad);
end

%%
fid = fopen('Readme.txt','w');
fprintf(fid,'Original folder name: %s \n',pwd);
fprintf(fid,'Date: %s \n\n\n',date);

fprintf(fid,'Model states and right hand side of the ODEs:\n');
for m=1:length(ar.model)
    for i=1:length(ar.model(m).xNames)
        fprintf(fid,'%s\t\t%s\n',ar.model(m).xNames{i},ar.model(m).fx{i});
    end
end
fprintf(fid,'\n\n');


fprintf(fid,'Observables and Errors:\n');
for m=1:length(ar.model)
    if isfield(ar.model(m),'data')
        for d=1:length(ar.model(m).data)
            for i=1:length(ar.model(m).data(d).fy)
                fprintf(fid,'%s:\t\tfy=%s\t\tfystd=%s\n',ar.model(m).data(d).yNames{i},ar.model(m).data(d).fy{i},ar.model(m).data(d).fystd{i});
            end
        end
    end
end
fprintf(fid,'\n\n');

fclose(fid);

diary Readme.txt
arPrint
diary off
movefile('Readme.txt',pfad);

fid = fopen([pfad,filesep,'Readme.txt'],'a');
fprintf(fid,'\n\n');
        
fprintf(fid,'Please write down a short documentation:\n');
fprintf(fid,'\n\n');
fclose(fid);

%% add path to git
nweg = length('arFramwork3/');
fid = fopen([pfad0(1:(end-nweg-1)),'.gitignore'],'a')
fprintf(fid,'%s\n',['arFramework3/Examples/_myOwnExamples/',name,'*']);
fclose(fid);

%%

edit([pfad,filesep,'Readme.txt']);


