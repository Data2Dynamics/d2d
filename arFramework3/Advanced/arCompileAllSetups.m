% setup_files = arCompileAllSetups(max_depth, mode)
%
% This function searches Setup files and executes the basic commands arInit,
% arLoadModel, arLoadData and arCompileAll to quickly compile serveral
% example models.
%
%       max_depth   If larger than 0, then the setup files be search
%                   recursively starting form the current directory.
%                   Default: Inf (recursively evaluates all setups found in
%                   deeper folders)
%       mode        'normal' (default: whole setup file)
%                   'core' only the most important commands, i.e. arInit,
%                   arLoadModel, arLoadData, arCompileAll
%
%      setup_files  array of found setup files
%
% One possibility to start it in the background (under linux):
% nohup matlab -nosplash < arCompileAllSetups_call.m > nohup.out &
%
%   Examples
%   arCompileAllSetups(2) % evaluates all Setups which are in pwd and 1 or
%                         % two folders deeper. This option might be
%                         % reasonable if executed in the folder
%                         % ....d2d/arFramework/Examples   

function setup_files = arCompileAllSetups(max_depth, mode)
    if ~exist('max_depth','var') || isempty(max_depth)
        max_depth = Inf;
    end
    if ~exist('mode','var') || isempty(mode)
        mode = 'normal';
    end
    
    testmode = false; % for code improvement
    
    if max_depth>0
        all_files = list_files_recursive(max_depth);
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
    
    pause_state = pause('query');
    pause('off');
    
    for i=1:length(setup_files)
        fprintf(fidlog,'%s will be executed ...\n\n' ,setup_files{i});
        [pathstr,name,ext] = fileparts(setup_files{i});
        % change path (if not there yet)
        if ~isempty(pathstr)
            cd(pathstr);
        end
        
        try
            switch mode
                case 'normal'
                    % evaluate the whole setup file
                    doEvalSetup(name);
                    if ~testmode
                        % ar is saved after executing the whole setup file
                        arSave(sprintf('arCompileAllSetups_%s', mode));
                    end
                case 'core'
                    % scan through the setup file and evaluate only the most important
                    % setup commands: arInit, arLoadModel, arLoadData, arCompileAll
                    global ar
                    did_compile = false;
                    fid = fopen([name,ext], 'r');
                    while (~feof(fid))
                        %             [str, fid] = arTextScan(fid, '%s', 'Delimiter', '', 'CommentStyle', '%');
                        str = fgets(fid);
                        %             disp(str);
                        
                        if ~isempty(str) && ischar(str)
                            str = strtrim(str);
                            if ~isempty(str) && ischar(str)
                                if ~isempty(regexp(str,'^arInit'))
                                    fprintf(fidlog,'%s\n',str);
                                    if ~testmode
                                        doEvalLine(str);
                                    end
                                elseif ~isempty(regexp(str,'^arLoadModel\('))
                                    fprintf(fidlog,'%s\n',str);
                                    if ~testmode
                                        doEvalLine(str);
                                    end
                                elseif ~isempty(regexp(str,'^arLoadData\('))
                                    fprintf(fidlog,'%s\n',str);
                                    if ~testmode
                                        doEvalLine(str);
                                    end
                                elseif ~isempty(regexp(str,'^arCompileAll'))
                                    fprintf(fidlog,'%s\n',str);
                                    did_compile = true;
                                    if ~testmode
                                        doEvalLine(str);
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
            pause(pause_state);
            warning('\n%s failed !!!!!!!!!', setup_files{i});
            warning(ERR.message);
            fprintf(fidlog,'\n%s failed !!!!!!!!! \n\n', setup_files{i});
            cd(pfad);
        end
    end
    
    fclose(fidlog);
    cd(pfad)
    pause(pause_state);
    
    end
    
    
    function doEvalSetup(fun_setup)
    % Evaluate a whole file (here: setup file) in its own workspace.
    % All local variables are protected from being overwritten or cleared.
    
    feval(fun_setup);
    end
    
    
    function doEvalLine(str)
    % Evaluate a string with limited scope (only affect ar struct).
    % All local variables are protected from being overwritten or cleared.
    
    global ar
    eval(str)
    end
    
    
    
    
    
    