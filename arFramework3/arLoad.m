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


if(isfield(ar,'isCompressed') && ar.isCompressed ~=0)
    arUncompress
end




function arUncompress(file)
% Inverse function to arCompress
% 
%   The sensitivities are filled again and the dimension infos are removed.
if(~exist('file','var'))
    file = '';
end

global ar
if(~isempty(file))
    load(file,'ar')
end

if(isfield(ar,'isCompressed') && ar.isCompressed ~=0)
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            ar.model(m).condition(c).suFineSimu = zeros(ar.model(m).condition(c).suFineSimu_dim);
            ar.model(m).condition(c).svFineSimu = zeros(ar.model(m).condition(c).svFineSimu_dim);
            ar.model(m).condition(c).sxFineSimu = zeros(ar.model(m).condition(c).sxFineSimu_dim);
            ar.model(m).condition(c).suExpSimu = zeros(ar.model(m).condition(c).suExpSimu_dim);
            ar.model(m).condition(c).svExpSimu = zeros(ar.model(m).condition(c).svExpSimu_dim);
            ar.model(m).condition(c).sxExpSimu = zeros(ar.model(m).condition(c).sxExpSimu_dim);
        end
        ar.model(m).condition = rmfield(ar.model(m).condition,'suFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'svFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'sxFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'suExpSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'svExpSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'sxExpSimu_dim');
    end
    
    try % for computers where simulation is not feasible
        arChi2(true);
    end
end
