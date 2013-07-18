% Stochastic Simulation Algorithm (SSA)
% see in Gillespie, The Journal of Physical Chemistry, 1977, 81(25)
%
% arSSA(m, scaling, nruns)
%
% m:        model index                             [1]
% scaling:  rescaling factor from species       
%           concentration to number of molecules    [ones(1,length(ar.model(m).x))]
% nruns:    number of runs                          [10]

function arSSA(m, scaling, nruns)

global ar

if(~exist('nruns','var'))
    nruns = 10;
end
ar.config.ssa_runs = nruns;
if(~exist('scaling','var') || isempty(scaling))
    scaling = ones(1,length(ar.model(m).x));
end

% update initial concentrations
arSimu(false, true);

% flux scaling
kscaling = zeros(1, length(ar.model(m).fv));
for j=1:length(ar.model(m).fv)
    qi = find(ar.model(m).N(:,j)<0);
    if(~isempty(qi))
        kscaling(j) = scaling(qi(1));
    else
        qi = find(ar.model(m).N(:,j)>0);
        if(~isempty(qi))
            kscaling(j) = scaling(qi(1));
        else
            error('reaction with empty N');
        end
    end
end

% setup arrays
for c=1:length(ar.model(m).condition)
    ar.model(m).condition(c).xFineSSA = nan(length(ar.model(m).condition(c).tFine), ...
        length(ar.model(m).x), nruns);
    ar.model(m).condition(c).xFineSSA_lb = nan(length(ar.model(m).condition(c).tFine), ...
        length(ar.model(m).x), nruns);
    ar.model(m).condition(c).xFineSSA_ub = nan(length(ar.model(m).condition(c).tFine), ...
        length(ar.model(m).x), nruns);
    
    if(isfield(ar.model(m).condition(c), 'tExp'))
        ar.model(m).condition(c).xExpSSA = nan(length(ar.model(m).condition(c).tExp), ...
            length(ar.model(m).x), nruns);
    end
    
    ar.model(m).condition(c).x0_ssa = ar.model(m).condition(c).xFineSimu(1,:);
    ar.model(m).condition(c).scale_x_ssa = scaling;
    ar.model(m).condition(c).scale_v_ssa = kscaling;
end
if(isfield(ar.model(m),'data'))
    for d=1:length(ar.model(m).data)
        ar.model(m).data(d).yFineSSA = nan(length(ar.model(m).data(d).tFine), ...
            length(ar.model(m).data(d).y), nruns);
        ar.model(m).data(d).yFineSSA_lb = nan(length(ar.model(m).data(d).tFine), ...
            length(ar.model(m).data(d).y), nruns);
        ar.model(m).data(d).yFineSSA_ub = nan(length(ar.model(m).data(d).tFine), ...
            length(ar.model(m).data(d).y), nruns);
        
        if(isfield(ar.model(m).data(d), 'tExp'))
            ar.model(m).data(d).yExpSSA = nan(length(ar.model(m).data(d).tExp), ...
                length(ar.model(m).data(d).y), nruns);
        end
    end
end

% call SSA mex function
fprintf('Gillespie simulation (%i runs)...', ar.config.ssa_runs); tic;
eval(sprintf('%s_SSA(ar);', ar.fkt));
fprintf('done (%g sec)\n', toc);

