% check systems setup
% addpath and configure sundials

function arCheck

if(verLessThan('symbolic', '5.5'))
	error('MUPAD symbolic toolbox version >= 5.5 required');
end

ar_path = strrep(which('arInit.m'),'/arInit.m','');
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

% % configure sundials 2.4.0
% if(exist([ar_path '/sundials-2.4.0'],'dir') == 0)
%     path_backup = cd;
%     cd(ar_path);
%     !tar -xvf sundials-2.4.0.tar
%     cd(path_backup);
% end
% if(exist([ar_path '/sundials-2.4.0/config.h'],'file') == 0)
%     path_backup = cd;
%     cd([ar_path '/sundials-2.4.0']);
%     !./configure
%     cd(path_backup);
% end

% configure sundials 2.5.0
if(exist([ar_path '/sundials-2.5.0'],'dir') == 0)
    path_backup = cd;
    cd(ar_path);
    !tar -xvf sundials-2.5.0.tar
    cd(path_backup);
end
if(exist([ar_path '/sundials-2.5.0/config.h'],'file') == 0)
    path_backup = cd;
    cd([ar_path '/sundials-2.5.0']);
    !./configure
    cd(path_backup);
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
	comment_string = input('\nPlease define a comment string (default is ''//'')\n-> ','s');
	if isempty(comment_string)
		comment_string = '//';
	end
	fprintf(fid,'\n%s%s%s','ar.config.comment_string = ''',comment_string,''';');
	fclose(fid);
	fprintf(1,'\n%s\n','arInitUser.m has been successfully created!');
    rehash path
end