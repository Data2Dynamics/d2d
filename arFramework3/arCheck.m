% check systems setup
% addpath and configure sundials

function docontinue = arCheck

if(verLessThan('symbolic', '5.5'))
	error('MUPAD symbolic toolbox version >= 5.5 required');
end

if(~ispc)
    ar_path = strrep(which('arInit.m'),'/arInit.m','');
else
    ar_path = strrep(which('arInit.m'),'\arInit.m','');
end

addpath(ar_path)

% load path of sub-directories
if(exist('pleInit','file') == 0)
    addpath([ar_path '/PLE2'])
end
if(exist('fileChooser','file') == 0)
    addpath([ar_path '/arTools'])
end
if(exist('JEInterface','file') == 0)
    addpath([ar_path '/EvA2/JEInterface'])
end

% uncompress and expand CVODES
if(exist([ar_path '/sundials-2.5.0'],'dir') == 0)
    if(~ispc)
        path_backup = cd;
        cd(ar_path);
        !tar -xvf sundials-2.5.0.tar
        cd(path_backup);
    else
        fprintf('Please uncompress and expand the CVODES sources\n%s\nand repeat.\n', [ar_path '\sundials-2.5.0.tar']);
        docontinue = false;
        return
    end
end

% write sundials_config.h
if(exist([ar_path '/sundials-2.5.0/include/sundials/sundials_config.h'],'file') == 0)
    fid = fopen([ar_path '/sundials-2.5.0/include/sundials/sundials_config.h'], 'W');
    fprintf(fid, '#define SUNDIALS_PACKAGE_VERSION "2.5.0"\n');
    fprintf(fid, '#define SUNDIALS_DOUBLE_PRECISION 1\n');
    fprintf(fid, '#define SUNDIALS_USE_GENERIC_MATH\n');
    fprintf(fid, '#define SUNDIALS_BLAS_LAPACK 0\n');
    fprintf(fid, '#define SUNDIALS_EXPORT\n');
    fclose(fid);
end

% EvA2 Toolbox
javaaddpath([ar_path '/EvA2/EvA2Base.jar'])

% check if arInitUser.m exists and create the file if necessary
if exist('arInitUser.m','file')==0
	fprintf(1,'\n\n%s\n%s\n\n','WARNING: arInitUser.m does not exist!','Creating the file...');
	user = '';
	while isempty(user)
		user = input('Please enter your full name (e.g. John Doe)\n-> ','s');
	end
	fid = fopen( [ar_path '/arInitUser.m'],'w');
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
