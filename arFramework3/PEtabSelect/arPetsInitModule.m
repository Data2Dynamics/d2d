global pets
petsInitMain();

function petsInitMain()

% global pets
global pets
pets = struct();

% start python and verify it is configured
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
