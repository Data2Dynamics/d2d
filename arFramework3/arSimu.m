% Simulate for current parameter settings
%
% arSimu(sensi, fine, dynamics)
%   sensi:          calculate sensitivities         [true]
%   fine:           fine grid for plotting          [false]
%   dynamics:       evaluate dynamics               [true]
% 
% or
%
% ar = arSimu(ar, sensi, fine, dynamics)
%   ar:             d2d model/data structure

function varargout = arSimu(varargin)

global ar

if(nargin==0 || ~isstruct(varargin{1}))
    qglobalar = true;
else
    ar = varargin{1};
    if(nargin>1)
        varargin = varargin(2:end);
    else
        varargin = {};
    end
    qglobalar = false;
end

if(~isempty(varargin) && ~isempty(varargin{1}))
    sensi = varargin{1};
else
    sensi = true;
end
if(length(varargin)>1 && ~isempty(varargin{2}))
    fine = varargin{2};
else
    fine = false;
end
if(length(varargin)>2 && ~isempty(varargin{3}))
    dynamics = varargin{3};
else
    dynamics = 0;
end

% If no sensitivities are requested, we always simulate (fast anyway)
if ( ~sensi )
    dynamics = 1;
end

% If dynamics are not requested, but sensitivities are, check whether we 
% already simulated these sensitivities.
% This code *only* prevents unnecessary simulation of sensitivities.
if ( ~dynamics && sensi )
    % Check cached config settings to see if they are still the same. If
    % not, then ar.cache.fine and ar.cache.exp get cleared.
    arCheckCache;
    
    % Check whether ar.cache.fine and exp are different from the ones we simulated last time
    % If not, we need to resimulate.
    if ( fine && ( ~isequal( ar.cache.fine(ar.qDynamic==1), ar.p(ar.qDynamic==1) ) ) )
        dynamics = 1;
    end
    if ( ~fine && ( ~isequal( ar.cache.exp(ar.qDynamic==1), ar.p(ar.qDynamic==1) ) ) )
        dynamics = 1;
    end
end

% We simulate sensitivities. Store the new 'last sensitivity simulation' parameters in the cache.
if ( dynamics && sensi )
    if ( fine )
        ar.cache.fine = ar.p;
    else
        ar.cache.exp = ar.p;
    end
end

% If we finally decide on not simulating the sensitivities, simulate only
% the model without sensitivities, since these could have been overwritten
% in the meantime
if ( ~dynamics )
    dynamics = 1;
    sensi = 0;
end

if(~isfield(ar,'p'))
    fprintf('ERROR: forgot arLink\n');
end
if(~isfield(ar.config,'useParallel'))
    ar.config.useParallel = true;
end
if(~isfield(ar.config,'fiterrors_correction'))
    ar.config.fiterrors_correction = 1;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

ar.stop = 0;

ss_presimulation = 0;
if ( isfield( ar, 'ss_conditions' ) && ( ar.ss_conditions ) )
    ss_presimulation = 1;
end

% propagate parameters
for m=1:length(ar.model)
    if ( ss_presimulation )
        for c=1:length(ar.model(m).ss_condition)       
            ar.model(m).ss_condition(c).status = 0;
            ar.model(m).ss_condition(c).pNum = ar.p(ar.model(m).ss_condition(c).pLink);
            ar.model(m).ss_condition(c).qLog10 = ar.qLog10(ar.model(m).ss_condition(c).pLink);
            ar.model(m).ss_condition(c).pNum(ar.model(m).ss_condition(c).qLog10 == 1) = ...
                10.^ar.model(m).ss_condition(c).pNum(ar.model(m).ss_condition(c).qLog10 == 1);
            ar.model(m).ss_condition(c).start = 0;
            ar.model(m).ss_condition(c).stop = 0;
            ar.model(m).ss_condition(c).stop_data = 0;
        end
    end
    
    for c=1:length(ar.model(m).condition)       
        ar.model(m).condition(c).status = 0;
        ar.model(m).condition(c).pNum = ar.p(ar.model(m).condition(c).pLink);
        ar.model(m).condition(c).qLog10 = ar.qLog10(ar.model(m).condition(c).pLink);
        ar.model(m).condition(c).pNum(ar.model(m).condition(c).qLog10 == 1) = ...
            10.^ar.model(m).condition(c).pNum(ar.model(m).condition(c).qLog10 == 1);
        ar.model(m).condition(c).start = 0;
        ar.model(m).condition(c).stop = 0;
        ar.model(m).condition(c).stop_data = 0;
    end
    
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            ar.model(m).data(d).pNum = ar.p(ar.model(m).data(d).pLink);
            ar.model(m).data(d).qLog10 = ar.qLog10(ar.model(m).data(d).pLink);
            ar.model(m).data(d).pNum(ar.model(m).data(d).qLog10 == 1) = ...
                10.^ar.model(m).data(d).pNum(ar.model(m).data(d).qLog10 == 1);
        end
    end
end

% initialize fine sensitivities
% this is very important, c code crashes otherwise!
if(fine && sensi)
    ar = initFineSensis(ar);
end

% do we have steady state presimulations?
if ( ss_presimulation )
    if ( sensi )
        ar = initSteadyStateSensis(ar);
    end
    feval(ar.fkt, ar, true, ar.config.useSensis && sensi, dynamics, false, 'ss_condition', 'ss_threads');
    
    % integration error ?
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).ss_condition)
            if(ar.model(m).ss_condition(c).status>0)
                error('arSimuCalc failed at %s for model %i, condition %i during pre-equilibration %i', ar.info.arsimucalc_flags{ar.model(m).ss_condition(c).status}, m, ar.model(m).ss_condition(c).src, c);
            elseif(ar.model(m).ss_condition(c).status<0)
                error('cvodes failed at %s for model %i, condition %i during pre-equilibration %i', ...
                    ar.info.cvodes_flags{abs(ar.model(m).ss_condition(c).status)}, m, ar.model(m).ss_condition(c).src, c);
            end
        end
    end
    
    for m = 1 : length( ar.model )
        % Map the steady states onto the respective target conditions
        for ssID = 1 : length( ar.model(m).ss_condition )
            targetConditions    = ar.model(m).ss_condition(ssID).ssLink;
            SSval               = ar.model(m).ss_condition(ssID).xFineSimu(end, :) + 0;     % + 0 is for R2013 compatibility
            if ( sensi )
                SSsens              = ar.model(m).ss_condition(ssID).sxFineSimu(end, :, :) + 0;
            end
            
            nStates = length(ar.model(m).x);
            
            % Which states do we change?
            ssStates = ar.model(m).ss_condition(ssID).ssStates;
            ssIgnore = ar.model(m).ss_condition(ssID).ssIgnore;
            
            % Copy the steady state values and sensitivities into the target
            % conditions taking into account any parameter order remapping
            for a = 1 : length( targetConditions )
                ar.model(m).condition(targetConditions(a)).modx_A(1,:) = 1 - ssStates;
                ar.model(m).condition(targetConditions(a)).modx_B(1,:) = SSval .* ssStates;
                if ( sensi )
                    ar.model(m).condition(targetConditions(a)).modsx_A(1,:,:) = ...
                        zeros(1,nStates,length(ar.model(m).condition(targetConditions(a)).p));
                    ar.model(m).condition(targetConditions(a)).modsx_B(1,:,:) = ...
                        SSsens(:,:,ar.model(m).condition(targetConditions(a)).ssParLink) + 0;
                    
                    % Certain sensitivities are not mappable from SS to target
                    % do not tamper with these.
                    if ~isempty( ar.model(m).condition(targetConditions(a)).ssUnmapped )
                        ar.model(m).condition(targetConditions(a)).modsx_A(1,:,ar.model(m).condition(targetConditions(a)).ssUnmapped) = 1;
                        ar.model(m).condition(targetConditions(a)).modsx_B(1,:,ar.model(m).condition(targetConditions(a)).ssUnmapped) = 0;
                    end
                    if ~isempty( ssIgnore )
                        ar.model(m).condition(targetConditions(a)).modsx_A(1,ssIgnore,:) = 1;
                        ar.model(m).condition(targetConditions(a)).modsx_B(1,ssIgnore,:) = 0;
                    end
                end
            end
        end
    end
end

% call mex function to simulate models
feval(ar.fkt, ar, fine, ar.config.useSensis && sensi, dynamics, false, 'condition', 'threads')

% integration error ?
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(ar.model(m).condition(c).status>0)
            error('arSimuCalc failed at %s', ar.info.arsimucalc_flags{ar.model(m).condition(c).status});
        elseif(ar.model(m).condition(c).status<0)
            error('cvodes failed at %s for model %i, condition %i. Trying arCheckSimu could be an option. ', ...
                ar.info.cvodes_flags{abs(ar.model(m).condition(c).status)}, m, c);
        end
    end
end

% manually transform sensitivities from normal to log10 for fine time points
if(fine && sensi && ar.config.useSensis)
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            for j=find(ar.qLog10(ar.model(m).condition(c).pLink)==1)
                ar.model(m).condition(c).suFineSimu(:,:,j) = ar.model(m).condition(c).suFineSimu(:,:,j) * ...
                    ar.model(m).condition(c).pNum(j) * log(10);
                ar.model(m).condition(c).sxFineSimu(:,:,j) = ar.model(m).condition(c).sxFineSimu(:,:,j) * ...
                    ar.model(m).condition(c).pNum(j) * log(10);
                ar.model(m).condition(c).szFineSimu(:,:,j) = ar.model(m).condition(c).szFineSimu(:,:,j) * ...
                    ar.model(m).condition(c).pNum(j) * log(10);
            end
        end
        if(isfield(ar.model(m), 'data'))
            for d=1:length(ar.model(m).data)
                for j=find(ar.qLog10(ar.model(m).data(d).pLink)==1)
                    ar.model(m).data(d).syFineSimu(:,:,j) = ar.model(m).data(d).syFineSimu(:,:,j) * ...
                        ar.model(m).data(d).pNum(j) * log(10);
                    ar.model(m).data(d).systdFineSimu(:,:,j) = ar.model(m).data(d).systdFineSimu(:,:,j) * ...
                        ar.model(m).data(d).pNum(j) * log(10);
                end
            end
        end
    end
end

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = {};
end

% (Re-)Initialize arrays for fine sensitivities with zeros
function ar = initFineSensis(ar)

for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = zeros(length(ar.model(m).data(d).tFine), ...
                length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
            ar.model(m).data(d).systdFineSimu = zeros(length(ar.model(m).data(d).tFine), ...
                length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
        end
    end
    for c = 1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu = zeros(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).u), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).svFineSimu = zeros(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).vs), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).sxFineSimu = zeros(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).x), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).szFineSimu = zeros(length(ar.model(m).condition(c).tFine), ...
            length(ar.model(m).z), length(ar.model(m).condition(c).p));
    end
end

function ar = initSteadyStateSensis(ar)

for m = 1:length(ar.model)
    for c = 1:length(ar.model(m).ss_condition)
        ar.model(m).ss_condition(c).suFineSimu = zeros(length(ar.model(m).ss_condition(c).tFine), ...
            length(ar.model(m).u), length(ar.model(m).ss_condition(c).p));
        ar.model(m).ss_condition(c).svFineSimu = zeros(length(ar.model(m).ss_condition(c).tFine), ...
            length(ar.model(m).vs), length(ar.model(m).ss_condition(c).p));
        ar.model(m).ss_condition(c).sxFineSimu = zeros(length(ar.model(m).ss_condition(c).tFine), ...
            length(ar.model(m).x), length(ar.model(m).ss_condition(c).p));       
        ar.model(m).ss_condition(c).szFineSimu = zeros(length(ar.model(m).ss_condition(c).tFine), ...
            length(ar.model(m).z), length(ar.model(m).ss_condition(c).p));
    end
end
