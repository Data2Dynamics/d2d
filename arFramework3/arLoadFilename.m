% load model struct and last ple results

function [ar, pleGlobals] = arLoadFilename(filename)

arCheck;

if(~exist('filename', 'var'))
    [~, filename] = fileChooser('./Results', 1, true);
end

S = load(['./Results/' filename '/workspace.mat']);
ar = S.ar;
fprintf('workspace loaded from file %s\n', filename);

pleGlobals = pleLoad(ar);
