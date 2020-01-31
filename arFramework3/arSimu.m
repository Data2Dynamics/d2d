%varargout = arSimu(varargin)
%
% arSimu simulates model dynamics and steady states for current parameter 
% settings
%
%   arSimu
%   arSimu([sensi], [fine], [dynamics])
%   ar = arSimu(ar, --)
% 
%   sensi       logical indicating whether to calculate sensitivities     [true]
%   fine        logical indicating whether to use fine grid for plotting  [false]
%   dynamics    logical indicating wheter to evaluate dynamics, if not 
%               simulate only steady states                               [true]
%   ar          (optional) d2d model/data structure
%
% Examples
%    arSimu(ar)    
%
% See also arPlot

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

% Computing scale parameters hierarchically requires the data grid, i.e. fine==false.
% So in the case of a fine grid, i.e. fine==true, we need a double pass:
% 1) simulate the model with the data grid to figure out the scale parameters and then
% 2) simulate the model with a fine grid, using the scale parameters resulting from the previous step.
% Here we carry out the 1st of these steps.
if isfield(ar.config,'useHierarchical') && ar.config.useHierarchical && fine
    arSimu(sensi,false,dynamics)
end

% If dynamics are not forced, check whether the dynamics of the last simulation
% were identical. If not, we have to resimulate.
if ( ~dynamics )
    % Check cached config settings to see if they are still the same. If
    % not, then cache storage gets cleared forcing resimulation.
    dynamics = arCheckCache;
    
    % Check whether dynamic parameters are different from the ones we 
    % simulated last time. If so, we need to resimulate!
    if ( ~dynamics )
        if ( fine )
            if ( ~isequal( ar.cache.fine(ar.qDynamic==1), ar.p(ar.qDynamic==1) ) || ( ar.cache.fineSensi ~= sensi ) )
                dynamics = 1;
            end
        else
            if ( ~isequal( ar.cache.exp(ar.qDynamic==1), ar.p(ar.qDynamic==1) ) || ( ar.cache.expSensi ~= sensi ) )
                dynamics = 1;
            end
        end
    end
end

% Make sure the required fields exist!
if ( isfield( ar.config, 'sensitivitySubset' ) && ( ar.config.sensitivitySubset == 1 ) )
    arSubsetSensi;
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

ss_presimulation = false;
if ( isfield( ar, 'ss_conditions' ) && ( ar.ss_conditions ) )
    ss_presimulation = true;
end

% propagate parameters
for m=1:length(ar.model)
    if ( ss_presimulation )
        for c=1:length(ar.model(m).ss_condition)   
            ar.model(m).ss_condition(c).x0_override = []; % remove initial condition overrides which may have been used for rootfinding
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
if(sensi)
    if ( fine )
        ar = initFineSensis(ar, dynamics);
    else
        ar = initExpSensis(ar, dynamics);
    end
end

% do we have steady state presimulations?
% options for steady state simulation:
%   1. Full simulation of everything (slow)
%       turboSSSensi = 0, rootFinding = 0
%   2. Fast simulation (simulate system without sensitivities, then determine sensis with implicit func theorem) (REQUIRES REDUCED SYSTEM)
%       turboSSSensi = 1, rootFinding = 0
%   3. Rootfinding; no guarantee that the determined steady state is the steady state the system would equilibrate to when multiple steady states exist (REQUIRES REDUCED SYSTEM)
%       rootFinding = 1
if ( ss_presimulation && dynamics )
    rootFinding = isfield( ar.config, 'rootfinding' ) && ( ar.config.rootfinding > 0.1 );
    if ( sensi || rootFinding )
        ar = initSteadyStateSensis(ar, dynamics);
    end
    if ( ~rootFinding )
        if ( isfield( ar.config, 'turboSSSensi' ) && ( ar.config.turboSSSensi == 1 ) )
            % Steady state determination by simulation without sensitivities and then determining them via implicit func theorem (only valid when conserved moieties have been removed from the model)
            for m=1:length(ar.model)
                fastSteadyState( m, sensi, dynamics );
            end
        else
            % Steady state determination by full simulation
            feval(ar.fkt, ar, true, ar.config.useSensis && sensi, dynamics, false, 'ss_condition', 'ss_threads', ar.config.skipSim);
        end
    else
        % Steady state determination by rootfinding
        for m=1:length(ar.model)
            for c=1:length(ar.model(m).ss_condition)
                [x, S] = arFindRoots(m, c, 'ss_condition', 1);
                ar.model(m).ss_condition(c).xFineSimu(end, :) = x + 0;
                if ( sensi )
                    ar.model(m).ss_condition(c).sxFineSimu(end, :, :) = S + 0;
                end
            end
        end
    end
    
    % integration error ?
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).ss_condition)
            if(ar.model(m).ss_condition(c).status>0)
                arCheckCache(1); % Invalidate cache so simulations do not get skipped
                nonEq = sprintf('%s ', ar.model(m).x{find(abs(ar.model(m).ss_condition(c).dxdt)>ar.config.eq_tol)});
                error('arSimuCalc failed at %s for model %i, condition %i during pre-equilibration %i.\nStates which failed to equilibrate:\n%s.\nCan attempt debugging using arPlotEquilibration.\n', ar.info.arsimucalc_flags{ar.model(m).ss_condition(c).status}, m, ar.model(m).ss_condition(c).src, c, nonEq);
            elseif(ar.model(m).ss_condition(c).status<0)
                arCheckCache(1); % Invalidate cache so simulations do not get skipped
                error('cvodes failed at %s for model %i, condition %i during pre-equilibration %i\nCan attempt debugging using arPlotEquilibration.\n', ...
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
                    
                    % Certain state sensitivities are not mappable from SS
                    % to target. If set, these should have an initial
                    % condition of zero for the sensitivities, since they
                    % do not exist in the reference state
                    if ~isempty( ar.model(m).condition(targetConditions(a)).ssUnmapped )
                        ar.model(m).condition(targetConditions(a)).modsx_A(1,:,ar.model(m).condition(targetConditions(a)).ssUnmapped) = 0;
                        ar.model(m).condition(targetConditions(a)).modsx_B(1,:,ar.model(m).condition(targetConditions(a)).ssUnmapped) = 0;
                    end
                    
                    % Sensitivities which are explicitly ignored are not
                    % affected
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
if ( isfield( ar.config, 'onlySS' ) && ( ar.config.onlySS == 1 ) )
    % Even if we only simulate steady states, we still need to propagate
    % the initial sensi to the observables.
    feval(ar.fkt, ar, fine, ar.config.useSensis && sensi, dynamics, false, 'condition', 'threads', 1)
else
    feval(ar.fkt, ar, fine, ar.config.useSensis && sensi, dynamics, false, 'condition', 'threads', ar.config.skipSim)
end

% integration error ?
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(ar.model(m).condition(c).status>0)
            arCheckCache(1); % Invalidate cache so simulations do not get skipped
            error('arSimuCalc failed at %s', ar.info.arsimucalc_flags{ar.model(m).condition(c).status});
        elseif(ar.model(m).condition(c).status<0)
            arCheckCache(1); % Invalidate cache so simulations do not get skipped
            error('cvodes failed at %s for model %i, condition %i. Trying arCheckSimu could be an option. ', ...
                ar.info.cvodes_flags{abs(ar.model(m).condition(c).status)}, m, c);
        end
    end
end

if(~fine)
    % arCalcRes_test;  % the test can be performed here or outside of arSimu
    % (see comments in arCalcRes_test.m)
    arCalcRes(sensi)

    %calculate y_scale_max
    if(ar.config.atolV)
        for im = 1:length(ar.model)
            for ic = 1:length(ar.model(im).condition)
                ar.model(im).condition(ic).atol(:) = 0;
                for id = ar.model(im).condition(ic).dLink
                    atol_tmp = squeeze(nanmax(nanmax(ar.model(im).data(id).y_scale,[],1),[],2));
                    if(~isequal(size(atol_tmp),size(ar.model(im).condition(ic).y_atol)))
                       warning('Problems in vector tolerances, did you load an old model?');                        
                    else
                        ar.model(im).condition(ic).y_atol = nanmax([ar.model(im).condition(ic).y_atol atol_tmp],[],2);
                    end
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

% We managed to complete the simulation
% Note that this has to happen at the end, since otherwise the simulation
% could have been terminated by the user in the meantime
if ( dynamics )
    if ( fine )
        ar.cache.fine               = ar.p + 0;
        ar.cache.fineSensi          = sensi + 0;
    else
        ar.cache.exp                = ar.p + 0;
        ar.cache.expSensi           = sensi + 0;
    end
end

% (Re-)Initialize arrays for fine sensitivities with zeros
function ar = initFineSensis(ar, dynamics)

for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = zeros(length(ar.model(m).data(d).tFine), ...
                length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
            ar.model(m).data(d).systdFineSimu = zeros(length(ar.model(m).data(d).tFine), ...
                length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
        end
    end
    if ( dynamics )
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
end

function ar = initExpSensis(ar, dynamics)

if ( dynamics )
    for m = 1:length(ar.model)
        if(isfield(ar.model(m), 'data'))
            for d = 1:length(ar.model(m).data)
                if ( ar.model(m).data(d).has_tExp )
                    ar.model(m).data(d).syExpSimu = zeros(length(ar.model(m).data(d).tExp), ...
                        length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
                    ar.model(m).data(d).systdExpSimu = zeros(length(ar.model(m).data(d).tExp), ...
                        length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
                end
            end
        end
    end
end

if ( isfield( ar.config, 'onlySS' ) && ( ar.config.onlySS == 1 ) )
	for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            if ( isfield( ar.model(m).condition(c), 'sxExpSimu' ) )
                ar.model(m).condition(c).sxExpSimu = zeros( size( ar.model(m).condition(c).sxExpSimu ) );
                ar.model(m).condition(c).svExpSimu = zeros( size( ar.model(m).condition(c).svExpSimu ) );
                ar.model(m).condition(c).szExpSimu = zeros( size( ar.model(m).condition(c).szExpSimu ) );
                ar.model(m).condition(c).suExpSimu = zeros( size( ar.model(m).condition(c).suExpSimu ) );
            end
        end
	end
end


function ar = initSteadyStateSensis(ar, dynamics)

if ( dynamics )
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
end

function fastSteadyState( m, sensi, dynamics )
    global ar;

    % Steady state determination by simulation
    %                 fine  sensi  dynamics  ssa    which condition field
    feval(ar.fkt, ar, true, false, dynamics, false, 'ss_condition', 'ss_threads', ar.config.skipSim);
	if ( ar.config.useSensis && sensi )
        for c = 1 : numel(ar.model(m).ss_condition)
            method = 1;
            dfdx = ar.model(m).ss_condition(c).dfdxNum + 0;
            dfdp = ar.model(m).ss_condition(c).dfdpNum + 0;
            
            if ( method == 2 )
                [Sx, r] = linsolve(-dfdx,dfdp); % For invertibility, model may not have conserved moieties
                        
                %if ( r < eps(1) )
                %    error( 'Model has conserved moieties or has not been sufficiently equilibrated. Fast equilibration result may be unreliable. Unless you know what you are doing, turn ar.config.turboSSSensi off by invoking ar.config.turboSSSensi = 0 or reduce the model prior to compilation (see help arReduce)' );
                %end
            else
                Sx = pinv(-dfdx)*dfdp;
                %Sx = dfdx.' \ dfdp;
            end
            ar.model(m).ss_condition(c).sxFineSimu(end,:,:) = Sx + 0;
            
            % Compute flux sensitivities
            dvdp = ar.model(m).ss_condition(c).dvdpNum;
            dvdx = ar.model(m).ss_condition(c).dvdxNum;
            Sv = dvdx * Sx + dvdp;
            ar.model(m).ss_condition(c).svFineSimu(end,:,:) = Sv + 0;
        end
	end