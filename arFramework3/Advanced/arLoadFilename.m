% load model struct and last ple results
%
% WARNING: 
% for some reason, after 
% > ar = arLoadFilename(filename)
% another
% > global ar
% has to be executed. 

function [arout, pleGlobals] = arLoadFilename(filename)

arCheck;

if(~exist('filename', 'var'))
    [~, filename] = fileChooser('./Results', 1, true);
end

S = load(['./Results/' filename '/workspace.mat']);
arout = S.ar;
fprintf('workspace loaded from file %s\n', filename);

pleGlobals = pleLoad(arout);
