% arCompileAllSetups
%
% setup_files = arCompileAllSetups(recursive)
%
%       recursive   should Setup files be search recursively starting form
%                   the current directory?
%                   Default: true
%
%       mode        'normal' (default: whole setup file)
%                   'core' only the most important commands, i.e. arInit,
%                   arLoadModel, arLoadData, arCompileAll
%
%   This function searches Setup files and executes the basic commands arInit,
%   arLoadModel, arLoadData and arCompileAll to quickly compile serveral
%   example models.
%
% One possibility to start it in the background (under linux):
% nohup matlab -nosplash < arCompileAllSetups_call.m > nohup.out &
%


function setup_files = arCompileAllSetups(recursive,mode)
if ~exist('recursive','var') || isempty(recursive)
    recursive = true;
end
if ~exist('mode','var') || isempty(mode)
    mode = 'normal';
end

testmode = false; % for code improvement

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

fidlog = fopen('arCompileAllSetups.log','w');

setup_files = all_files(bol);
fprintf('The following setup files were found and will be subsequently used for compiling: \n\n');
fprintf(fidlog,'The following setup files were found and will be subsequently used for compiling: \n\n');
for i=1:length(setup_files)
    fprintf('%s\n',setup_files{i});
    fprintf(fidlog,'%s\n',setup_files{i});
end
fprintf('\n\n');
fprintf(fidlog,'\n\n');

% create parallel pool if not yet existing:
p = gcp('nocreate');
if isempty(p) && ~testmode
    parpool('local')
end

pfad = pwd;

for i=1:length(setup_files)
    fprintf(fidlog,'%s will be executed ...\n\n' ,setup_files{i});
    [pathstr,name,ext] = fileparts(setup_files{i});
    cd(pathstr);
    
    try
        switch mode
            case 'normal' % every line
                eval(name);
                if ~testmode
                    eval(sprintf('arSave(''arCompileAllSetups_%s'');',mode));
                end
            case 'core' % only the most improtant setup commands
                fid = fopen([name,ext], 'r');
                while (~feof(fid))
                    %             [str, fid] = arTextScan(fid, '%s', 'Delimiter', '', 'CommentStyle', '%');
                    str = fgets(fid);
                    %             disp(str);
                    
                    if ~isempty(str) && ischar(str)
                        str = strtrim(str);
                        did_compile = false;
                        if ~isempty(str) && ischar(str)
                            if ~isempty(regexp(str,'^arInit'))
                                fprintf(fidlog,'%s\n',str);
                                if ~testmode
                                    eval(str);
                                end
                            elseif ~isempty(regexp(str,'^arLoadModel\('))
                                fprintf(fidlog,'%s\n',str);
                                if ~testmode
                                    eval(str);
                                end
                            elseif ~isempty(regexp(str,'^arLoadData\('))
                                fprintf(fidlog,'%s\n',str);
                                if ~testmode
                                    eval(str);
                                end
                            elseif ~isempty(regexp(str,'^arCompileAll'))
                                fprintf(fidlog,'%s\n',str);
                                did_compile = true;
                                if ~testmode
                                    eval(str);
                                end
                            else
                                
                            end
                        end
                    end
                end % lines
                
                if did_compile
                    if ~testmode
                        eval(sprintf('arSave(''arCompileAllSetups_%s'');',mode));
                    end
                end
                
                fprintf(fidlog,'\n\n');
                fclose(fid);
            otherwise
                error('mode %s is not implemented.',model)
        end  % switch mode
    catch ERR
        warning('\n%s failed !!!!!!!!! \n\n' ,setup_files{i});
        disp(lasterror);
        fprintf(fidlog,'\n%s failed !!!!!!!!! \n\n' ,setup_files{i});
        cd(pfad);
    end
end

fclose(fidlog);
cd(pfad)






