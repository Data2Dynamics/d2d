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
            ar.model(m).condition(c).uFineSimu = zeros(ar.model(m).condition(c).uFineSimu_dim);
            ar.model(m).condition(c).vFineSimu = zeros(ar.model(m).condition(c).vFineSimu_dim);
            ar.model(m).condition(c).xFineSimu = zeros(ar.model(m).condition(c).xFineSimu_dim);
            ar.model(m).condition(c).zFineSimu = zeros(ar.model(m).condition(c).zFineSimu_dim);
            
            ar.model(m).condition(c).suFineSimu = zeros(ar.model(m).condition(c).suFineSimu_dim);
            ar.model(m).condition(c).svFineSimu = zeros(ar.model(m).condition(c).svFineSimu_dim);
            ar.model(m).condition(c).sxFineSimu = zeros(ar.model(m).condition(c).sxFineSimu_dim);
            ar.model(m).condition(c).szFineSimu = zeros(ar.model(m).condition(c).szFineSimu_dim);
        end
        ar.model(m).condition = rmfield(ar.model(m).condition,'uFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'vFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'xFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'zFineSimu_dim');
        
        ar.model(m).condition = rmfield(ar.model(m).condition,'suFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'svFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'sxFineSimu_dim');
        ar.model(m).condition = rmfield(ar.model(m).condition,'szFineSimu_dim');
    end
    
    try % for computers where simulation is not feasible
        arChi2(true);
    end
end
