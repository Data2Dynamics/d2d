% load model struct and last ple results

arCheck;

global ar
global pleGlobals

[~, filename] = fileChooser('./Results', 1, true);

Stmpload = load(['./Results/' filename '/workspace.mat']);
ar = Stmpload.ar;

fprintf('workspace loaded from file %s\n', filename);

pleGlobals = pleLoad(ar);

clear Stmpload
clear filename