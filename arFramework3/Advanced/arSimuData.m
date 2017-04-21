% Simulate data for current parameter settings
%
% arSimuData(m, jplot, tpoints, randomseed)
%   tpoints:    time points for simulation        
%   m:          model index                    
%   jplot:      plot index
%   randomseed  random seed for noise generation

function arSimuData(m, jplot, tpoints, randomseed)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end 

if ~exist('tpoints','var')
    tpoints = [];
end

if ~exist('randomseed','var') || isempty(randomseed)
    randomseed = 'shuffle';
end

% set random seed
if(exist('rng','file')~=0)
    if(exist('randomseed','var') && ~isempty(randomseed))
        ar.simudata_seed = randomseed;
        rng(randomseed);
    else
        rng(randomseed);
        rngsettings = rng;
        ar.simudata_seed = rngsettings.Seed;
    end
end

if(~exist('jplot','var') || isempty(jplot))
    jplot = [];
end

if(~exist('m','var') || isempty(m))
    for jm=1:length(ar.model)
        arSimuData(jm, jplot, tpoints, randomseed);
    end
    return
else % m index provided:
    if(isempty(jplot))
        for jplot=1:length(ar.model(m).plot)
            arSimuData(m, jplot, tpoints, randomseed);
        end
        return
    elseif(length(jplot)>1)
        for jjplot=jplot
            arSimuData(m, jjplot, tpoints, randomseed);
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

% remember existing parameters
tmp_p = ar.p;
tmp_lb = ar.lb;
tmp_ub = ar.ub;
tmp_qFit = ar.qFit;
tmp_qLog10 = ar.qLog10;

arLink(true);

ar.p = tmp_p;
ar.lb = tmp_lb;
ar.ub = tmp_ub;
ar.qFit = tmp_qFit;
ar.qLog10 = tmp_qLog10;


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
% ar.model(m).data(d).yExpStd = zeros(length(tpoints),
% length(ar.model(m).data(d).y)); % don't overwrite exp. errors
ar.model(m).data(d).yExpSimu = zeros(length(tpoints), length(ar.model(m).data(d).y));
ar.model(m).data(d).ystdExpSimu = zeros(length(tpoints), length(ar.model(m).data(d).y));

function simulate_data(m, d)
global ar

if all(size(ar.model(m).data(d).yExp) == size(ar.model(m).data(d).yExpStd))
    yExpStdSimu = ar.model(m).data(d).yExpStd;
    nosd = isnan(ar.model(m).data(d).yExpStd); % don't overwrite exp. errors by the error model
    yExpStdSimu(nosd) = ar.model(m).data(d).ystdExpSimu(nosd);
else
    yExpStdSimu = ar.model(m).data(d).ystdExpSimu;
end

ar.model(m).data(d).yExp = ar.model(m).data(d).yExpSimu + ...
    randn(size(ar.model(m).data(d).yExpSimu)) .* yExpStdSimu;

ar.model(m).data(d).tLim = ar.model(m).data(d).tLimExp;
ar.model(m).data(d).tLim(2) = round(max(ar.model(m).data(d).tExp)*1.1);