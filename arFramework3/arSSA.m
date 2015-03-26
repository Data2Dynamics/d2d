% Stochastic Simulation Algorithm (SSA)
% see in Gillespie, The Journal of Physical Chemistry, 1977, 81(25)
%
% arSSA(nruns,scaling)
%
% nruns:    number of runs                          [10]
% scaling:  rescaling factor from species           [{ones}]
%           concentration to number of molecules
    
function arSSA(nruns, scaling)

global ar

if(~exist('nruns','var'))
    nruns = 10;
end
ar.config.ssa_runs = nruns;

% update initial concentrations
arSimu(false, true);

for m=1:length(ar.model)
    if(~exist('scaling','var'))
        scaling{m} = ones(1,length(ar.model(m).x));
    end


    % flux scaling
    kscaling = zeros(1, length(ar.model(m).fv));
    for j=1:length(ar.model(m).fv)
        qi = find(ar.model(m).N(:,j)<0);
        if(~isempty(qi))
            kscaling(j) = scaling{m}(qi(1));
        else
            qi = find(ar.model(m).N(:,j)>0);
            if(~isempty(qi))
                kscaling(j) = scaling{m}(qi(1));
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
        
        ar.model(m).condition(c).zFineSSA = nan(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).z), nruns);
        ar.model(m).condition(c).zFineSSA_lb = nan(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).z), nruns);
        ar.model(m).condition(c).zFineSSA_ub = nan(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).z), nruns);
        
        if(isfield(ar.model(m).condition(c), 'tExp'))
            ar.model(m).condition(c).xExpSSA = nan(length(ar.model(m).condition(c).tExp), ...
                length(ar.model(m).x), nruns);
            ar.model(m).condition(c).zExpSSA = nan(length(ar.model(m).condition(c).tExp), ...
                length(ar.model(m).z), nruns);
        end
        
        ar.model(m).condition(c).x0_ssa = ar.model(m).condition(c).xFineSimu(1,:);
        ar.model(m).condition(c).scale_x_ssa = scaling{m};
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
end

% fprintf('Stochastic simulation (%i runs)...', ar.config.ssa_runs); tic;

% call mex function to simulate models
feval(ar.fkt, ar, true, false, true, true, 'condition', 'threads')

% fprintf('done (%g sec)\n', toc);

