% check systems setup

function arCheck

if(verLessThan('symbolic', '5.5'))
	error('MUPAD symbolic toolbox version >= 5.5 required');
end

ar_path = strrep(which('arInit.m'),'/arInit.m','');

% load path of sub-directories
if(exist('ple','file') == 0)
    addpath([ar_path '/PLE2'])
end
if(exist('fileChooser','file') == 0)
    addpath([ar_path '/arTools'])
end
if(exist('JEInterface','file') == 0)
    addpath([ar_path '/EvA2/JEInterface'])
end

% configure sundials
if(exist([ar_path '/sundials-2.4.0'],'dir') == 0)
    path_backup = cd;
    cd(ar_path);
    !tar -xvf sundials-2.4.0.tar
    cd(path_backup);
end
if(exist([ar_path '/sundials-2.4.0/config.h'],'file') == 0)
    path_backup = cd;
    cd([ar_path '/sundials-2.4.0']);
    !./configure
    cd(path_backup);
end

% EvA2 Toolbox
javaaddpath([ar_path '/EvA2/EvA2Base.jar'])