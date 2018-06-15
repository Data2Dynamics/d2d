function checkstr = arChecksumFitting(arStruct)
if nargin==0
    global ar
    arStruct = ar;
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




% This function defines all fields in arStruct to be checked for testing whether
% fitting properties are equal. Model properties are not covered! This
% should be covered by arStruct.fkt (i.e. by the global checksum) and by
% arCheckSumParameter
function a = initCheckFields(arStruct)

a = struct;

% fields to be checkt in arStruct.config
fn = intersect(fieldnames(arStruct.config),{'add_c','atol','atolV','atolV_Sens','eq_rtol','eq_step_factor','eq_tol', ...
    'fastEquilibration','fiterrors','init_eq_step','max_eq_steps','maxsteps',...
    'maxstepsize','maxtol','nCVRestart','nfine_dr_method','optimizer','rootfinding',...
    'rtol','sensiSkip','sensitivitySubset','skipSim','ssa_min_tau','ssa_runs',...
    'steady_state_constraint','turboSSSensi','turboSplines','useEvents','useFitErrorCorrection',...
    'useFitErrorMatrix','useJacobian','useMS','useSensiRHS','useSensis','useSparseJac',...
    'useTolSwitching','useTolTrustPar','optimizerStep'});

for i=1:length(fn)
    a.(fn{i}) = arStruct.config.(fn{i});
end


a.config.optim = struct;
% fields to be checkt in arStruct.config.optim
fnoptim = setdiff(fieldnames(arStruct.config.optim),{'Display','DerivativeCheck','Diagnostics','FunValCheck','OutputFcn','PlotFcns','UseParallel'});
for i=1:length(fnoptim)
    a.config.optim.(fnoptim{i}) = arStruct.config.optim.(fnoptim{i});
end

