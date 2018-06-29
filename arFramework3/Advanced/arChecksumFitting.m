%   arStruct        if instead of the global ar, the checksum should be
%                   evaluated for another struct, then it is provided as
%                   first argument
% 
%  saveEvaluatedFields  Default: false
%                   if true, then a workspace is saved in folder Checksums
%                   containing the field which are evaluated for
%                   calculationg the checksum

function checkstr = arChecksumFitting(arStruct, saveEvaluatedFields)
global ar
if ~exist('arStruct','var') || isempty(arStruct)
    arStruct = ar;
end
if ~exist('saveEvaluatedFields','var') || isempty(saveEvaluatedFields)
    saveEvaluatedFields = false;
end

if saveEvaluatedFields
    arCopy = struct;
end

if ~isfield(arStruct,'config')
    checkstr = '';
    return
elseif ~isfield(arStruct.config,'optim')
    checkstr = '';
    return    
end

a = initCheckFields(arStruct);
checksum = arAddToCheckSum(a);

h = typecast(checksum.digest,'uint8');
checkstr = dec2hex(h)';
checkstr = checkstr(:)';

clear checksum


if saveEvaluatedFields
    arSaveChecksumCopy(a,'fitting',checkstr);
end



% This function defines all fields in arStruct to be checked for testing whether
% fitting properties are equal. Model properties are not covered! This
% should be covered by arStruct.fkt (i.e. by the global checksum) and by
% arCheckSumParameter
function a = initCheckFields(arStruct)

a = struct;
a.config = struct;

% fields to be checkt in arStruct.config
fn = intersect(fieldnames(arStruct.config),{'add_c','atol','atolV','atolV_Sens','eq_rtol','eq_step_factor','eq_tol', ...
    'fastEquilibration','fiterrors','init_eq_step','max_eq_steps','maxsteps',...
    'maxstepsize','maxtol','nCVRestart','nfine_dr_method','optimizer','rootfinding',...
    'rtol','sensiSkip','sensitivitySubset','skipSim','ssa_min_tau','ssa_runs',...
    'steady_state_constraint','turboSSSensi','turboSplines','useEvents','useFitErrorCorrection',...
    'useFitErrorMatrix','useJacobian','useMS','useSensiRHS','useSensis','useSparseJac',...
    'useTolSwitching','useTolTrustPar','optimizerStep','user_residual_fun'});

for i=1:length(fn)
    a.config.(fn{i}) = arStruct.config.(fn{i});
end


a.config.optim = struct;
% fields to be checkt in arStruct.config.optim
fnoptim = setdiff(fieldnames(arStruct.config.optim),{'Display','DerivativeCheck','Diagnostics','FunValCheck','OutputFcn','PlotFcns','UseParallel'});
for i=1:length(fnoptim)
    a.config.optim.(fnoptim{i}) = arStruct.config.optim.(fnoptim{i});
end


