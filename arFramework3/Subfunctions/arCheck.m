% docontinue = arCheck
% arCheck checks the system's setup, sets Matlab paths via addpath and
% configure sundials
%
%   docontinue      true if the function is finished, false otherwise

function docontinue = arCheck
docontinue = false;

symbtool = ver('symbolic');
if(~isempty(symbtool) && verLessThan('symbolic', '5.5'))
    error('MUPAD symbolic toolbox version >= 5.5 required');
end

ar_path = fileparts(which('arInit.m'));

% if the path to d2d-folder has been specified with a local path, then turn it to
% the global path:
arFolderName = strsplit(ar_path,filesep);
matpath = path;
matpath = strsplit(matpath,pathsep);
if ~any(strcmp(ar_path,matpath))
    warning('It seems that you specified the path to d2d using local paths, e.g. by addpath(''arFramework3''). The path is now switch to a gloabl one.')
    path_id = ~cellfun(@isempty,regexp(matpath,[arFolderName{end} '$']));
    if(~all(path_id)==0)
        rmpath(matpath{path_id})
    end
    addpath(ar_path)
end

% add all subfolders of arFramework3 folder to MATLAB search path
% load path of sub-directories
if(exist('arSIAInit','file') == 0)
    addpath([ar_path '/IA'])
end
if(exist('arPLEInit','file') == 0)
    addpath([ar_path '/PLE'])
end
if(exist('InitVPL','file') == 0)
    addpath([ar_path '/VPL'])
end
if(exist('Init2DPL','file') == 0)
    addpath([ar_path '/PL2D'])
    addpath([ar_path '/PL2D/Subfunctions_autofix'])
    addpath([ar_path '/PL2D/InterX'])
end
if(exist('arHelpStruct','file') == 0)
    addpath([ar_path '/Help'])
end
if(exist('doPPL','file') == 0)
    addpath([ar_path '/PPL'])
end
if(exist('fileChooser','file') == 0)
    addpath([ar_path '/MatlabTools'])
end
if(exist('l1Init','file') == 0)
    addpath([ar_path '/L1'])
end
if(exist('arChi2s','file') == 0)
    addpath([ar_path '/Advanced'])
end
if(exist('arMC3','file') == 0)
    addpath([ar_path '/MCMC'])
end
if(exist('arSimuCalc.c','file') == 0)
    addpath([ar_path '/Ccode'])
end
if(exist('arChi2Cluster','file') == 0)
    addpath([ar_path '/ClusterFunctions'])
end
if(exist('arParseModel','file') == 0)
    addpath([ar_path '/Deprecated'])
end
if(exist('arEvaluate','file') == 0)
    addpath([ar_path '/Development'])
end
if(exist('arToPython','file') == 0)
    addpath([ar_path '/d2d-presenter'])
end
if(exist('chemist.sty','file') == 0)
    addpath([ar_path '/Latex'])
end
if(exist('model_template.def','file') == 0)
    addpath([ar_path '/ProjectTemplate'])
end
if(exist('arFitLhsBwCluster','file') == 0)
    addpath([ar_path '/BwGrid']);
end
if(exist('arNEB','file') == 0)
    addpath([ar_path '/NEB']);
end
if(exist('arExportPEtab','file') == 0)
    addpath([ar_path '/ImportExport']);
end
if(exist('model_template_HillFunctions.def','file') == 0)
    addpath([ar_path '/ProjectTemplate']);
end
if(exist('arPetsInitModule','file') == 0)
    addpath([ar_path '/PEtabSelect']);
end


warning('off','MATLAB:rmpath:DirNotFound')

% removes Examples folder and subfolders of arFramework3 from the MATLAB
% search path to avoid loading data from those examples for accidentially
% identical file names
ids_Examples = strfindcell(matpath,[ar_path '/Examples']);
rmpath(strjoin(matpath(ids_Examples),':'));

% removes PM folder and subfolders of arFramework3
if(exist([ar_path,'/PM'], 'dir'))
    ids_PM = strfindcell(matpath,[ar_path '/PM']);
    rmpath(strjoin(matpath(ids_PM),':'));
    arFprintf(2, 'arCheck.m: rm PM from matlab path\n');
end

warning('on','MATLAB:rmpath:DirNotFound')


[has_git, is_repo] = arCheckGit(ar_path);

% check if submodules have been pulled from github
submodules = {'matlab2tikz','https://github.com/matlab2tikz/matlab2tikz/zipball/3a1ee10';...
    'parfor_progress','https://github.com/dumbmatter/parfor_progress/zipball/fcb9d86';...
    'plot2svg','https://github.com/jschwizer99/plot2svg/zipball/839a919';...
    'export_fig','https://github.com/altmany/export_fig/zipball/d8131e4';...
    'Ceres/ceres-solver','https://github.com/ceres-solver/ceres-solver/zipball/8ea86e1';...
    'NL2SOL','https://github.com/JoepVanlier/mexNL2SOL/zipball/daabaac';...
    };
for jsm = 1:length(submodules)
    submodule = submodules{jsm,1};
    submod_dir = [ar_path '/ThirdParty/' submodule];
    url = submodules{jsm,2};
    if( (exist(submod_dir,'dir')==7 && isempty(ls(submod_dir))) || (~exist(submod_dir,'file')) )
        arFprintf(2,'Downloading submodule "%s" from github...',submodule);
        if(has_git && is_repo)
            % fetch submodule via git
            library_path = getenv('LD_LIBRARY_PATH');
            setenv('LD_LIBRARY_PATH', '');
            old_path = pwd;
            cd(ar_path);
            if(isunix)
                system('git submodule update --init --recursive >/dev/null 2>&1');
            else
                system('git submodule update --init --recursive >nul 2>&1');
            end
            cd(old_path);
            setenv('LD_LIBRARY_PATH', library_path);
        else
            % fetch submodule as zip file
            if(exist(submod_dir,'dir'))
                rmdir(submod_dir);
            end
            fname = [submod_dir '.zip'];
            arDownload(url, fname);
            dirnames = unzip(fname, [ar_path filesep 'ThirdParty']);
            dirnames = unique(cellfun(@fileparts,dirnames,'UniformOutput',false));
            dirname = dirnames{1};
            movefile(dirname, submod_dir);
            delete(fname);
        end
        arFprintf(2,' done!\n');
    end
end


% path of third party software
if(exist('JEInterface','file') == 0)
    addpath([ar_path '/ThirdParty/EvA2/JEInterface'])
    javaaddpath([ar_path '/ThirdParty/EvA2/EvA2Base.jar'])
end
if(exist('suptitle','file') == 0)
    addpath([ar_path '/ThirdParty/BlandAltman'])
end
if(exist('export_fig','file') == 0)
    addpath([ar_path '/ThirdParty/export_fig'])
end
if(exist('plot2svg','file') == 0)
    if exist([ar_path '/ThirdParty/plot2svg/src'],'dir')
        addpath([ar_path '/ThirdParty/plot2svg/src'])
    end
end
if(exist('matlab2tikz','file') == 0)
    if exist([ar_path '/ThirdParty/matlab2tikz/src'],'dir')
        addpath([ar_path '/ThirdParty/matlab2tikz/src'])
    end
end
if(exist('parfor_progress','file') == 0)
    addpath([ar_path '/ThirdParty/parfor_progress'])
end
if (exist('compileNL2SOL', 'file') == 0)
    addpath([ar_path '/ThirdParty/NL2SOL'])
end
if (exist('mota', 'file') == 0)
    addpath([ar_path '/ThirdParty/MOTA'])
end
if (exist('arTRESNEI', 'file') == 0)
    addpath([ar_path '/ThirdParty/TRESNEI'])
end
if (exist('compileCeres', 'file') == 0)
    addpath([ar_path '/ThirdParty/Ceres'])
end
if (exist('TranslateSBML', 'file') == 0)
    addpath([ar_path '/ThirdParty/libSBML'])
end
if (exist('fminsearchbnd', 'file') == 0)
    addpath([ar_path '/ThirdParty/FMINSEARCHBND'])
end
if (exist('STRSCNE', 'file') == 0)
    addpath([ar_path '/ThirdParty'])
end
if (exist('YAMLMatlab', 'file') == 0)
    addpath([ar_path '/ThirdParty/YAMLMatlab'])
end
if (exist('StrucID', 'file') == 0)
    addpath([ar_path '/ThirdParty/StrucID'])
end
if (exist('StrucID/adimat', 'file') == 0)
    addpath([ar_path '/ThirdParty/StrucID/adimat'])
end
if (exist('strike-goldd', 'file') == 0)
    addpath([ar_path '/ThirdParty/strike-goldd'])
end

savepath

%% CVODES

% uncompress and expand CVODES
if(exist([ar_path '/ThirdParty/sundials-2.6.1'],'dir') == 0)
    path_backup = cd;
    cd([ar_path '/ThirdParty']);
    untar('sundials-2.6.1.tar');
    cd(path_backup);
end

% write sundials_config.h
if(exist([ar_path '/ThirdParty/sundials-2.6.1/include/sundials/sundials_config.h'],'file') == 0)
    fid = fopen([ar_path '/ThirdParty/sundials-2.6.1/include/sundials/sundials_config.h'], 'W');
    if(fid==-1)
        error('could not write file %s!', [ar_path '/ThirdParty/sundials-2.6.1/include/sundials/sundials_config.h']),
    end
    fprintf(fid, '#define SUNDIALS_PACKAGE_VERSION "2.6.1"\n');
    fprintf(fid, '#define SUNDIALS_DOUBLE_PRECISION 1\n');
    fprintf(fid, '#define SUNDIALS_USE_GENERIC_MATH\n');
    fprintf(fid, '#define SUNDIALS_BLAS_LAPACK 0\n');
    fprintf(fid, '#define SUNDIALS_EXPORT\n');
    fclose(fid);
end

%% SuiteSparse 4.2.1

% uncompress and expand KLU solver
if(exist([ar_path '/ThirdParty/KLU-1.2.1'],'dir') == 0)
    path_backup = cd;
    cd([ar_path '/ThirdParty']);
    untar('KLU-1.2.1.tar');
    cd(path_backup);
end


%% check Windows libraries for pthread-win32
if(ispc)
    %     if(exist(['.\pthreadGC2_',mexext,'.dll'],'file')==0)
    try
        copyfile([ar_path '\ThirdParty\pthreads-w32_2.9.1\dll\' mexext '\pthreadGC2.dll'], 'pthreadGC2.dll');
    catch ERR  % occurs (and can be ignored), if dll has been copied previously, is still loaded and therefore replacement is blocked by Windows OS
        disp(ERR.message)
    end
    %     end
    %     if(exist(['.\pthreadVC2_',mexext,'.dll'],'file')==0)
    try
        copyfile([ar_path '\ThirdParty\pthreads-w32_2.9.1\dll\' mexext '\pthreadVC2.dll'], 'pthreadVC2.dll');
    catch ERR  % occurs (and can be ignored), if dll has been copied previously, is still loaded and therefore replacement is blocked by Windows OS
        disp(ERR.message)
    end
    %     end
end

%% user name

% check if arInitUser.m exists and create the file if necessary
if exist('arInitUser.m','file')==0
    fprintf(1,'\n%s\n\n','Welcome to Data 2 Dynamics Software');
    user = '';
    while isempty(user)
        user = input('Please enter your full name (e.g. John Doe)\n-> ','s');
    end
    arCreateInitUser(user);
end

docontinue = true;
