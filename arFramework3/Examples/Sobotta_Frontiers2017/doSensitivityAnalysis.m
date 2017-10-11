addpath SensitivityAnalysis

if ( exist( 'BigFiles/Final_Profiles/Final_Core_profiles.mat', 'file' ) )
    % Backward compatibility
    pleGlobals      = load('BigFiles/Final_Profiles/Final_Core_profiles');
    pleGlobals      = pleGlobals.ple;
    
    AveragedSensitivityAnalysis;
    disp( 'The result was stored in /figures/LPSA' );
else
    SensitivityAnalysis;
end