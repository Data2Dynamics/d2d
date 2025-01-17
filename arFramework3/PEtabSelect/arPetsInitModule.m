% ARPETSINITMODULE Initialize the model selection module "PEtab Select"
%
% If the PEtab Select package is installed and loaded correctly, a global
% struct `pets` is created. It contains the PEtab Select module.
%
% If there are any issues with the installation, an error is thrown.
% For installation instructions for the PEtab Select module with d2d, see:
% https://github.com/Data2Dynamics/d2d/wiki/Model-selection-with-PEtab-Select

global pets
petsInitMain();

function petsInitMain()

% global pets
global pets
pets = struct();

% initialize pets.history
pets.history = createEmptyHistory();

% start python and verify it is configured correctly
pe = pyenv;
if isempty(pe.Version)
    error('Python is not configured in MATLAB. Use the pyenv function to set up Python.');
end

verifyPythonVersion(pe, 3, 10);

% try to load petab_select package
try 
    pets.module = py.importlib.import_module('petab_select');
    pets.plot = py.importlib.import_module('petab_select.plot');
    fprintf('The petab_select package is installed and loaded correctly.\n');
catch
    error('The petab_select package is not installed. Please install it using pip.');
end

end


function verifyPythonVersion(pe, expectMajor, expectMinor)
% VERIFYPYTHONENVIRONEMENT check version requirements of python environment
%
% INPUTS:
%   pe (Python environment object)
%   expectMajor (double): Expected major version
%   expectMinor (double): Expected minor version
%
% Throws an error if the Python version is less than the expected version.


% Get the Python version
pythonVersion = pe.Version;

% Split the version string into major, minor, and patch
versionParts = split(pythonVersion, '.');
if length(versionParts) < 2
    error('Could not parse Python version: %s', pythonVersion);
end

% Convert to numeric
major = str2double(versionParts{1});
minor = str2double(versionParts{2});

% Check version
if major > expectMajor || (major == expectMajor && minor >= expectMinor)
    fprintf('Python version is %i.%i or higher.\n', expectMajor, expectMinor);
else
    error('Python version is less than %i.%i. Current version: %s', ...
        expectMajor, expectMinor, pythonVersion);
end

end


function history = createEmptyHistory()
% create raw structure for history of model selection

history = struct();

% struct for currently useful models
history.current = struct();
history.current.predecessor = '';
history.current.iteration_best = '';
history.current.allTime_best = '';

% struct for calibrated models
% subfields are initialized at first call of 'addModel2History'
history.calibrated = struct();

% lists for model subspace yaml files and names
history.PetabYamlFiles = {};
history.PetabNames = {};

end