% check systems setup
% addpath and configure sundials

function docontinue = arCheck

symbtool = ver('symbolic');
if(~isempty(symbtool) && verLessThan('symbolic', '5.5'))
	error('MUPAD symbolic toolbox version >= 5.5 required');
end

if(~ispc)
    ar_path = strrep(which('arInit.m'),'/arInit.m','');
else
    ar_path = strrep(which('arInit.m'),'\arInit.m','');
end

warning('off','MATLAB:rmpath:DirNotFound')

% add all subfolders of arFramework3 folder to MATLAB search path
rmpath(genpath( [ar_path,'/PM']))
disp('arCheck.m: rm PM from matlab path');


% removes Examples folder and subfolders of arFramework3 from the MATLAB
% serach path to avoid loading data from those examples for accidentially
% identical file names
rmpath(genpath([ar_path '/Examples']))

warning('on','MATLAB:rmpath:DirNotFound')

% load path of sub-directories
if(exist('pleInit','file') == 0)
    addpath([ar_path '/PLE2'])
end
if(exist('doPPL','file') == 0)
    addpath([ar_path '/PPL'])
end
if(exist('fileChooser','file') == 0)
    addpath([ar_path '/arTools'])
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
    addpath([ar_path '/ThirdParty/export_fig-master'])
end
if(exist('plot2svg','file') == 0)
    addpath([ar_path '/ThirdParty/plot2svg-20141130/src'])
end
if(exist('matlab2tikz','file') == 0)
    addpath([ar_path '/ThirdParty/matlab2tikz-matlab2tikz-722609f/src'])
end
if(exist('parfor_progress','file') == 0)
    addpath([ar_path '/ThirdParty/parfor_progress'])
end
    
%% CVODES

% uncompress and expand CVODES
if(exist([ar_path '/sundials-2.6.1'],'dir') == 0)
    path_backup = cd;
    cd(ar_path);
    untar('sundials-2.6.1.tar');
    cd(path_backup);
end

% write sundials_config.h
% if(exist([ar_path '/sundials-2.5.0/include/sundials/sundials_config.h'],'file') == 0)
%     fid = fopen([ar_path '/sundials-2.5.0/include/sundials/sundials_config.h'], 'W');
if(exist([ar_path '/sundials-2.6.1/include/sundials/sundials_config.h'],'file') == 0)
    fid = fopen([ar_path '/sundials-2.6.1/include/sundials/sundials_config.h'], 'W');
    if(fid==-1)
        error('could not write file %s!', [ar_path '/sundials-2.6.1/include/sundials/sundials_config.h']),
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
if(exist([ar_path '/KLU-1.2.1'],'dir') == 0)
    path_backup = cd;
    cd(ar_path);
    untar('KLU-1.2.1.tar');
    cd(path_backup);
end


%% check Windows libraries for pthread-win32
if(ispc)
    if(exist('C:\Windows\pthreadGC2.dll','file')==0 || exist('C:\Windows\pthreadVC2.dll','file')==0)
        if(exist('.\pthreadGC2.dll','file')==0)
            copyfile([ar_path '\pthreads-w32_2.9.1\dll\' mexext '\pthreadGC2.dll'], 'pthreadGC2.dll');
        end
        if(exist('.\pthreadVC2.dll','file')==0)
            copyfile([ar_path '\pthreads-w32_2.9.1\dll\' mexext '\pthreadVC2.dll'], 'pthreadVC2.dll');
        end
    end
end

%% user name

% check if arInitUser.m exists and create the file if necessary
if exist('arInitUser.m','file')==0
	fprintf(1,'\n\n%s\n%s\n\n','WARNING: arInitUser.m does not exist!','Creating the file...');
	user = '';
	while isempty(user)
		user = input('Please enter your full name (e.g. John Doe)\n-> ','s');
	end
	fid = fopen([ar_path '/arInitUser.m'],'w');
    if(fid==-1)
        error('could not write file %s!', [ar_path '/arInitUser.m']),
    end
	fprintf(fid,'%s\n','% initialize user settings');
	fprintf(fid,'\n%s\n','function arInitUser');
	fprintf(fid,'\n%s\n','global ar');
	fprintf(fid,'\n%s%s%s','ar.config.username = ''',user,''';');
	fprintf(fid,'\n%s%s%s','ar.config.comment_string = ''//'';');
	fclose(fid);
	fprintf(1,'\n%s\n','arInitUser.m has been successfully created!');
    rehash path
end

docontinue = true;
