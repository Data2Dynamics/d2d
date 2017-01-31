function flag = arChangeVersion(version)

old_path = pwd;
ar_path = fileparts(which('arInit'));
cd(ar_path)

if(~exist('version','var') || strcmp(version,'current'))
    flag = system(sprintf('git checkout -f master'))==0;
else
    flag = system(sprintf('git checkout -f %s',version))==0;
end
cd(old_path);
