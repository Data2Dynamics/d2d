% Simulate data for current parameter settings
%
% arSimuData(m, jplot, tpoints)
%   tpoints:    time points for simulation        
%   m:          model index                    
%   jplot:      plot index                    

function arSimuData(m, jplot, tpoints)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if ~exist('tpoints','var')
    tpoints = [];
end

if(~exist('jplot','var') || isempty(jplot))
    jplot = [];
end

if(~exist('m','var') || isempty(m))
    for jm=1:length(ar.model)
        arSimuData(jm, jplot, tpoints);
    end
    return
else % m index provided:
    if(isempty(jplot))
        for jplot=1:length(ar.model(m).plot)
            arSimuData(m, jplot, tpoints);
        end
        return
    elseif(length(jplot)>1)
        for jjplot=jplot
            arSimuData(m, jjplot, tpoints);
        end     
        return
    end
end

ds = ar.model(m).plot(jplot).dLink;
if(isempty(tpoints))
    tpoints = ar.model(m).data(ds(1)).tExp;
end

% set time point and clear arrays
for d=ds
    assigne_new_timepoints(tpoints, m, d);
end
arLink(true);

% simulate model
arSimu(false,false);
ar.pTrue = ar.p;

% simulate data
for d=ds
    simulate_data(m, d);
end
arLink(true);


function assigne_new_timepoints(tpoints, m, d)
global ar
ar.model(m).data(d).tExp = sort(tpoints(:));
ar.model(m).data(d).yExp = zeros(length(tpoints), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpStd = zeros(length(tpoints), length(ar.model(m).data(d).y));
ar.model(m).data(d).yExpSimu = zeros(length(tpoints), length(ar.model(m).data(d).y));
ar.model(m).data(d).ystdExpSimu = zeros(length(tpoints), length(ar.model(m).data(d).y));

function simulate_data(m, d)
global ar
ar.model(m).data(d).yExp = ar.model(m).data(d).yExpSimu + ...
    randn(size(ar.model(m).data(d).yExpSimu)) .* ar.model(m).data(d).ystdExpSimu;
ar.model(m).data(d).yExpStd = ar.model(m).data(d).ystdExpSimu;
ar.model(m).data(d).tLim = ar.model(m).data(d).tLimExp;
