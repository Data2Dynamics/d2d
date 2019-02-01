% arResample([m], [d])
% 
% Resample data in ar.model(m).data(d) for current parameter settings
% 
%   m      index of model to resample data.              [all models]
%   d      index of data to resample data.               [all data]                   

function arResample(m, d)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if(~exist('m','var'))
    for jm=1:length(ar.model)
        for jd=1:length(ar.model(jm).data)
            arResample(jm, jd)
        end
    end
    return
end
if(~exist('d','var'))
    for jd=1:length(ar.model(m).data)
        arResample(m, jd)
    end
    return
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

% simulate model
arSimu(false,false);
ar.pTrue = ar.p;

% simulate data
qnan = isnan(ar.model(m).data(d).yExp);
if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors ~= -1) || ...
        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(m,d)~=-1) )
    ar.model(m).data(d).yExp = ar.model(m).data(d).yExpSimu + ...
        randn(size(ar.model(m).data(d).yExpSimu)) .* ar.model(m).data(d).ystdExpSimu;
else
    ar.model(m).data(d).yExp = ar.model(m).data(d).yExpSimu + ...
        randn(size(ar.model(m).data(d).yExpSimu)) .* ar.model(m).data(d).yExpStd;
end
ar.model(m).data(d).yExp(qnan(:)) = nan;
