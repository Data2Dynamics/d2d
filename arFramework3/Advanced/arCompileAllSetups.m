% nohup matlab -nosplash < arCompileAllSetups_call.m > nohup.out &


function setup_files = arCompileAllSetups(recursive)
if ~exist('recursive','var') || isempty(recursive)
    recursive = true;
end


if recursive
    all_files = list_files_recursive;
else
    d = dir;
    files = {d.name};
    files = files(find(~[ d.isdir]));
    all_files = setdiff(files,{'.','..'});
end

if ~iscell(all_files)
    all_files = {all_files};
end

bol = false(size(all_files));
for i=1:length(all_files)
    [pathstr,name,ext] = fileparts(all_files{i});
    bol(i) = ~isempty(regexp(name,'^[sS]etup*')) && strcmp(ext,'.m')==1;
end

setup_files = all_files(bol);
% cdisp(setup_files)


try
   delete(gcp('nocreate'))
end
try
   parpool('local')
end

pfad = pwd;
fidlog = fopen('arCompileAllSetups.log','w');

for i=1:length(setup_files)
    fprintf(fidlog,'%s will be executed ...\n\n' ,setup_files{i});
    [pathstr,name,ext] = fileparts(setup_files{i});
    cd(pathstr);
    
    try
        fid = fopen([name,ext], 'r');
        while (~feof(fid))
            [str, fid] = arTextScan(fid, '%s\n', 1, 'CommentStyle', '%');
            
            if ~isempty(str) && iscell(str)                
                str = strtrim(str{1});
                if ~isempty(str) && ischar(str{1})
                    if sum(~cellfun(@isempty,regexp(str,'arInit')))>0
                        fprintf(fidlog,'%s\n',str{1});
                        %         eval(str);
                    elseif sum(~cellfun(@isempty,regexp(str,'arLoadModel\(')))>0
                        fprintf(fidlog,'%s\n',str{1});
                        %         eval(str);
                    elseif sum(~cellfun(@isempty,regexp(str,'arLoadData\(')))>0
                        fprintf(fidlog,'%s\n',str{1});
                        %         eval(str);
                    elseif sum(~cellfun(@isempty,regexp(str,'arCompileAll')))>0
                        fprintf(fidlog,'%s\n',str{1});
                        %         eval(str);
                    else
                        
                    end
                else
%                     fprintf('%s\n',str{1});
                end
            end
        end
        fprintf(fidlog,'\n\n');
        fclose(fid);
    catch
        warning('\n%s failed !!!!!!!!! \n\n' ,setup_files{i});
        disp(lasterror);
        fprintf(fidlog,'\n%s failed !!!!!!!!! \n\n' ,setup_files{i});
        cd(pfad);
    end
end

fclose(fidlog);
cd(pfad)






