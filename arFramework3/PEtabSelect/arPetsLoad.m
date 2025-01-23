function arPetsLoad(name, PetsOutputDir)
% ARPETSLOAD loads the PEtab Select struct "pets" from file.
%
% The python module petab_select cannot be saved to/loaded from a .mat file
% Therefore, only the relevant fields are saved and loaded.
%
% INPUT:
%   name: Name of the file to load the struct from.
%   PetsOutputDir: Directory to load the file from.
%
% SEE: arPetsSave

arguments
    name char = 'petsSave';
    PetsOutputDir char = fullfile(pwd(), 'PEtabSelect', 'Results');
end

% initialize new pets struct
arPetsInitModule;

% load the struct from file
loadPath = fullfile(PetsOutputDir, name);
load(loadPath, "petsSave");

% assign the filed to the global variable
pets.history = petsSave.history;
pets.problem = petsSave.problem;

end
