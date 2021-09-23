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

command = [
sprintf(append('petab_select candidates \\ \n ',  ...
'-y model_selection\', filesep, 'selection_problem.yaml \\ \n ', ...
'-s output\', filesep, 'state.dill \\  \n ',...
'-o output\', filesep, 'models.yaml  \\ \n ', ...
'-m brute_force \\ \n ', ...
'-l 3') )
] ;
[status,cmdout] = system(command)


%% Import models.yaml written by petab_select script


%% Run d2d fits & calculate criterion values with d2d functions




%% Write calibrated first iteration .yaml including criterion values

%% use petab_select candidates -b calibrated first iteration.yaml with method (forward selection...)

%% Final output 

end