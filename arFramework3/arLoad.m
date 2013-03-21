% load model struct and last ple results

function arLoad(filename)

arCheck;

if(~exist('filename', 'var'))
    [~, filename] = fileChooser('./Results', 1, true);
end

load(['./Results/' filename '/workspace.mat']);
fprintf('workspace loaded from file %s\n', filename);
