% This function can be used to automatically update .gitignore if the set
% of examples has changed.
% The function automatically detects all folder in arFramework3/Examples/
function write_gitignore

fprintf('New .gitignore file will be written in the working directory ...\n\n')

fid = fopen('.gitignore','w');

% ignore cases independent on examples:
lines = {'# This file has been generated using arFramework3/arTools/write_gitignore.m'
    '#######################################'
    ''
    'syntax: glob'
    'arFramework3/arInitUser.m'
    'arFramework3/sundials-2.4.0/'
    'arFramework3/sundials-2.5.0/'
    'arFramework3/sundials-2.6.1/'
    'arFramework3/KLU-1.2.1/'
    '.DS_Store'
    ''
    '*.dat'
    ''
    'arFramework3/TestSuite'
    ''
    'arFramework3/sbml-test-cases-2014-10-22'
    'arFramework3/sbml-test-cases-2014-10-22.zip'
    ''
    'arFramework3/l1/trdog/trdog.m'
    ''};


% ignore cases for each example:
example_folder = strrep(which('arInit'),'arInit.m','Examples');
d = dir(example_folder);
folders = {d.name};
folders = folders([d.isdir]);
folders = setdiff(folders,{'.','..'});

for i=1:length(folders)
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/Compiled/'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/SBML/'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/arSimuCalc*'];
    lines{end+1} = ['arFramework3/Examples/',folders{i},'/arCheckParallel*'];
    lines{end+1} = '';
end

for i=1:length(lines)
    fprintf('%s\n',lines{i});
    fprintf(fid,'%s\n',lines{i});
end

fclose(fid);

fprintf('\n\n .gitignore has been written. Please check the file by hand.\n');
fprintf('If approved, replace the old file in the folder arFramework3 by hand.\n\n')