% Simulate data for current parameter settings
%
% arSimuData(tpoints, m, d)
%   tpoints:    time points for simulation        
%   m:          model index                    
%   d:          data index                    

function arSimuData(tpoints, m, d)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if(~exist('m','var'))
    for jm=1:length(ar.model) 
        arSimuData(tpoints, jm);
    end
    return
else
    if(~exist('d','var'))
        for jd=1:length(ar.model(m).data)
            arSimuData(tpoints, m, jd);
        end
        return
    end
end

ar.model(m).data(d).tExp = sort(tpoints(:));
ar.model(m).data(d).yExp = zeros(length(tpoints), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = zeros(length(tpoints), length(ar.model(m).data(d).y));

arLink(true);

% simulate model
ar.model(m).data(d).yExpSimu(:) = 0;
ar.model(m).data(d).ystdExpSimu(:) = 0;
arSimu(false,false);
ar.pTrue = ar.p;

% simulate data
ar.model(m).data(d).yExp = ar.model(m).data(d).yExpSimu + ...
    randn(size(ar.model(m).data(d).yExpSimu)) .* ar.model(m).data(d).ystdExpSimu;

ar.model(m).data(d).tLim = ar.model(m).data(d).tLimExp;
arLink(true);
