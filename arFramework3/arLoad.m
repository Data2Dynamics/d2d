% load model struct and last ple results
% 
%   filename    either the name of the result folder or its number

function arLoad(filename)
arCheck;

global ar
global pleGlobals

% set the two variables also as global in the command line workspace:
evalin('base','clear ar pleGlobals');  
evalin('base','global ar pleGlobals');  

if(~exist('filename', 'var') || isempty(filename))
    [~, filename] = fileChooser('./Results', 1, true);
elseif(isnumeric(filename)) % filename is the file-number
    [~, ~, file_list] = fileChooser('./Results', 1, -1);    
    filename = file_list{filename};
end

Stmpload = load(['./Results/' filename '/workspace.mat']);
ar = Stmpload.ar;

fprintf('workspace loaded from file %s\n', filename);

try
    pleGlobals = pleLoad(ar);
catch
    fprintf(1,'No valid PLE workspace found!\n');
    clear pleGlobals;
end




