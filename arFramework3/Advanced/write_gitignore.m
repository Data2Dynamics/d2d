% This function can be used to automatically update .gitignore if the set
% of examples has changed.
% The function automatically detects all folder in arFramework3/Examples/
function write_gitignore

fprintf('New .gitignore file will be written in the working directory ...\n\n')
ar_path = fileparts(which('arInit.m'));
nr_slash = strfind(ar_path,'/');
ar_path = ar_path(1:nr_slash(end));
fid = fopen([ar_path '.gitignore'],'w');

% ignore cases independent on examples:
lines = {'# This file has been generated using arFramework3/Advanced/write_gitignore.m'
    '#######################################'
    ''
    'syntax: glob'
    'arFramework3/arInitUser.m'
    'arFramework3/sundials-2.4.0/'
    'arFramework3/sundials-2.5.0/'
    'arFramework3/sundials-2.6.1/'
    'arFramework3/ThirdParty/sundials-2.6.1/'
    'arFramework3/KLU-1.2.1/'
    'arFramework3/ThirdParty/KLU-1.2.1/'
    'arFramework3/ThirdParty/Ceres/ceresd2d.mex*'
    '.DS_Store'
    ''
    '*.dat'
    ''
    'arFramework3/l1/trdog/trdog.m'
    ''
    'arFramework3/Examples/Dream6_L1/Data/simu*'
    'arFramework3/Examples/Dream6_L1/Results/20*'
    ''
    'arFramework3/Examples/_myOwnExamples/'
    ''};


% ignore cases for each example:
example_folder = strrep(which('arInit'),'arInit.m','Examples');
d = dir(example_folder);
subfolders = {'Biomodels','ToyModels'};

folders = {d.name};
folders = folders([d.isdir]);
folders = setdiff(folders,{'.','..'});
folders = setdiff(folders,subfolders);

for i=1:length(folders)
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/Compiled/'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/SBML/'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/arSimuCalc*'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/arCheckParallel*'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/Results/20*'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/pthreadVC2.dll'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/pthreadGC2.dll'];
    lines{end+1} = '';
end

for j = 1:length(subfolders)
    d = dir([example_folder '/' subfolders{j}]);
    folders = {d.name};
    folders = folders([d.isdir]);
    folders = setdiff(folders,{'.','..'});

    for i=1:length(folders)
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/Compiled/'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/SBML/'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/arSimuCalc*'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/arCheckParallel*'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/Results/20*'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/pthreadVC2.dll'];
        lines{end+1} = ['arFramework3/Examples/',subfolders{j},'/',folders{i},'/pthreadGC2.dll'];
        lines{end+1} = '';
    end
end

for i=1:length(lines)
    fprintf('%s\n',lines{i});
    fprintf(fid,'%s\n',lines{i});
end

fclose(fid);

fprintf('\n\n .gitignore has been written in d2d core folder. Please check the file by hand!\n');