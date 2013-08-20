% load model struct and last ple results

arCheck;

global ar
global pleGlobals

[~, filename] = fileChooser('./Results', 1, true);

S = load(['./Results/' filename '/workspace.mat']);
ar = S.ar;

fprintf('workspace loaded from file %s\n', filename);

pleGlobals = pleLoad(ar);

clear S
clear filename