function arPetsSave(name, PetsOutputDir)
% ARPETSSAVE saves the PEtab Select struct "pets" to file.
%
% The python module petab_select cannot be saved to/loaded from a .mat file
% Therefore, only the relevant fields are saved and loaded.
%
% INPUT:
%   name: Name of the file to save the struct to.
%   PetsOutputDir: Directory to save the file to.
%
% SEE: arPetsLoad

arguments
    name char = 'petsSave';
    PetsOutputDir char = fullfile(pwd(), 'PEtabSelect', 'Results');
end

% access the global pets struct
global pets

% create a struct with the relevant fields
petsSave = struct();
petsSave.history = pets.history;
petsSave.problem = pets.problem;

% save the struct to file
savePath = fullfile(PetsOutputDir, name);
save(savePath, "petsSave")

end