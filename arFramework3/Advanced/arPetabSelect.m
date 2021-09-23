%% Work in Progress

function arPetabSelect()
%% Parse function input


%% Write selection problem yaml with input


%% Check if petab_select is installed in system commands and give error if not yet installed

command = ['petab_select --help'];
[status,~] = system(command);
if status ~= 0
    error(sprintf('petab_select could not be called as a system command.\nPlease make sure to install petab_select using pip install petab_select on your machine.'))
end

%% Run petab_select code with system commands

% Create output folder? Currently still missing
if ~isfolder('output')
    mkdir('output')
end
command = [
append('petab_select candidates ',  ...
' -y selection_problem.yaml ', ...
' -s output', filesep, 'state.dill',...
' -o output', filesep, 'models.yaml', ...
' -m brute_force', ...
' -l 3')
] ;
[status,cmdout] = system(command)


%% Import models.yaml written by petab_select script


%% Run d2d fits & calculate criterion values with d2d functions




%% Write calibrated first iteration .yaml including criterion values

%% use petab_select candidates -b calibrated first iteration.yaml with method (forward selection...)

%% Final output 

end