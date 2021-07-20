global ar

strikegoldd_path = fileparts(mfilename('fullpath'));
addpath(strikegoldd_path)
addpath(fullfile(strikegoldd_path,'functions'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_functions'))
addpath(fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_models'))


ar.ia.addedpaths{1} = strikegoldd_path;
ar.ia.addedpaths{2} = fullfile(strikegoldd_path,'functions');
ar.ia.addedpaths{3} = fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_functions');
ar.ia.addedpaths{4} = fullfile(strikegoldd_path,'functions',filesep,'aux_Lie_symmetry',filesep,'lie_sym_models');

