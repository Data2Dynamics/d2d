% [ar2] = arCalcMerit([sensi], [pTrial], [dynamics], [doSimu])
% [ar2] = arCalcMerit([ar], [sensi], [pTrial], [dynamics], [doSimu])
%
% This function updates the fields representing objective functions (and
% parts of it) used for fitting for the currently chosen parameters.
%
%   sensi          propagate sensitivities         [false]
%                  this argument is passed to arSimu
%   pTrial         trial parameter of fitting
%   dynamics       force evaluation of dynamics    [false]
%   doSimu         should arSimu be called         [true]
%   ar             d2d model/data structure
%
%   ar2            d2d model/data structure with updated objective functions
%
% This function replaces arChi2 and is called by arFit. The fact that the
% function accepts serval calls is because of historical developments. It
% is not beautiful but required to not break compatibility.
%
% This function calls
%   - arSimu which in turn calls arCalcRes
%   - arCollectRes
%
% Possible calls:
% >> arCalcMerit        then:
%               qglobalar = true
%               sensi = true
%
% >> arCalcMerit(ar)    here, the global ar is overwritten by the argument
%               qglobalar = false
%               sensi = true
%
% >> arCalcMerit(ar,sensi)
% >> arCalcMerit(sensi)
%               like arCalcMerit, but sensi can be set
%
% >> arCalcMerit(ar,sensi,ptrial)
% >> arCalcMerit(sensi,ptrial)
%               like arCalcMerit(ar,sensi), but
%               silent = true
%               ar.p is set to ptrial
%               this is one possiblity to set 'silent' to true. The other
%               one is:
% >> arCalcMerit(sensi,[])
%
% >> arCalcMerit(sensi,ptrial,dynamics)
% >> arCalcMerit(ar,sensi,ptrial,dynamics)
%               dynamics is passed to arSimu, otherwise arSimu is called
%               without this argument, i.e. with its default
%
% >> arCalcMerit(sensi,ptrial,dynamics,doSimu)
% >> arCalcMerit(ar,sensi,ptrial,dynamics,doSimu)
%           doSimu  can be set to false, then the residuals are calculated
%           without updating the model trajectories (e.g. if ar.qFit or ar.
%           model.data.qFit) has been changed.
%
% Example (call used by arFit)
% ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
%
% See also arGetMerit, arChi2


function varargout = arCalcMerit(varargin)

global ar %#ok<*GVMIS>

nargs = nargin;
% The possiblity providing ar as an argument and to use of qglobalar==0 is
% obsolete because the gloal "ar" is overwritten anyway in arSimu
% Implementation due to backwards compability:
if nargs>0 && isstruct(varargin{1})
    ar = varargin{1};
    varargin = varargin(2:end);  % changing the meaning of varargin is not nicely implemented
    nargs = nargs-1;
    qglobalar = false;
else
    qglobalar = true;
end

if nargs>=1 && ~isempty(varargin{1})
    sensi = varargin{1} && ar.config.useSensis;
else
    sensi = true && ar.config.useSensis;
end

if nargs>=2 && ~isempty(varargin{2})
    pTrial = varargin{2};
    ar.p(ar.qFit==1) = pTrial + 0;
    silent = true;
elseif nargs==2 && isempty(varargin{2})
    silent = true;
else
    silent = false;
end

if nargs>=3 && ~isempty(varargin{3})
    dynamics = varargin{3};
else
    dynamics = [];
end

if nargs>=4 && ~isempty(varargin{4})
    doSimu = varargin{4};
else
    doSimu = true;
end

if(~isfield(ar, 'fevals'))
    ar.fevals = 0;
end
ar.fevals = ar.fevals + 1;


atol = ar.config.atol;
rtol = ar.config.rtol;
maxsteps = ar.config.maxsteps;
qPositiveX = {ar.model(:).qPositiveX};
allowNegativeX = 0;

if(~isfield(ar.config, 'nCVRestart') || isnan(ar.config.nCVRestart))
    nCVRestart = 10;
else
    nCVRestart = ar.config.nCVRestart;
end

for i = 1:nCVRestart
    
    if ~doSimu
        has_error = false;
        break
    end
    
    try
        if(qglobalar)  % since ar is overwritten anyway in arSimu, the possiblity to use of qglobalar obsolete
            arSimu(sensi, false, dynamics);
        else
            ar = arSimu(ar, sensi, false, dynamics);
        end
        
        has_error = false;
        break
        
    catch error_id
        
        has_error = true;
        if(~silent)
            arFprintf(1, 'Integration error.\n');
            disp(error_id.message);
        end
        
        % in the following cases we have no fix and exit the loop:
        % undefined function
        if strcmp(error_id.identifier,'MATLAB:UndefinedFunction')
            break
        end
        % no typical CVODES error flag -> there probably was another error
        % (there are also the CVODES flags 1, 2, and 99
        %  but currently they are not handled specifically)
        for m=1:length(ar.model)
            if all([ar.model(m).condition(:).status] >= 0)
                break
            end
        end
        
        % heuristic methods to try another time
        if i < nCVRestart
            
            increaseMaxsteps = 0;
            increaseTols = 0;
            
            for m=1:length(ar.model)
                
                % CV_TOO_MUCH_WORK
                if any([ar.model(m).condition(:).status] == -1)
                    
                    if any(ar.model(m).qPositiveX(:))
                        % some states have to positive?
                        % -> temporarily allow negative states
                        ar.model(m).qPositiveX(:) = 0;
                        allowNegativeX = 1;
                        
                    else
                        % no positivity restraints?
                        % -> increase max number of steps
                        increaseMaxsteps = 1;
                    end
                    
                end
                
                % other CVODES error, e.g., -2 = CV_TOO_MUCH_ACC (to much accuracy)
                if any([ar.model(m).condition(:).status] < -1)
                    increaseTols = 1;
                end
                
            end
            
            % console output
            arFprintf(1, 'New Attempt (%d/%d).\n', i+1, nCVRestart)
            if allowNegativeX
                arFprintf(1, 'Allow negative states.\n')
            end
            if increaseMaxsteps
                % why 20%?
                ar.config.maxsteps = (1.2)*ar.config.maxsteps;
                arFprintf(1, 'Increase maxsteps by 20%% (now: %0.2e).\n', ar.config.maxsteps)
            end
            if increaseTols
                % why 5%?
                ar.config.atol = (1.05)*ar.config.atol;
                ar.config.rtol = (1.05)*ar.config.rtol;
                arFprintf(1, 'Increase tolerances by 5%% (now: atol=%0.2e, rtol=%0.2e).\n', ...
                    ar.config.atol, ar.config.rtol)
            end
            
        end
    end
end

% reset parameters
ar.config.atol = atol;
ar.config.rtol = rtol;
ar.config.maxsteps = maxsteps;

for m=1:length(ar.model)
    ar.model(m).qPositiveX = qPositiveX{m};
    for c=1:length(ar.model(m).condition)
        if(isfield(ar.model(m).condition(c), 'xExpSimu'))
            if(sum((min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) ./ arRange(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),1) < -ar.config.rtol) & (min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) < -ar.config.atol)) > 0)
                arFprintf(1, 'Negative state in model %d condition %d detected that is defined as positive! Double-check model definition!\nPlot negative states by calling ar.model(%d).qPositiveX(:) = 0; with subsequent arPlot call.\n',m,c,m)
            end
        end
    end
end

arCollectRes(sensi);

% set Inf for errors
if(has_error)
    ar.res(:) = Inf;
    ar.chi2 = Inf;
    ar.chi2err = Inf;
    ar.chi2fit = Inf;
end

% calculate first order optimality criterion
if(sensi)
    res = [ar.res ar.constr];
    sres = [];
    if(~isempty(ar.sres))
        sres = ar.sres(:, ar.qFit==1);
    end
    if(~isempty(ar.sconstr))
        sres = [sres; ar.sconstr(:, ar.qFit==1)];
    end
    g = -2*res*sres; % gradient
    if(~isempty(g))
        onbound = [my_equals(ar.p(ar.qFit==1),ar.ub(ar.qFit==1)); my_equals(ar.p(ar.qFit==1),ar.lb(ar.qFit==1))];
        exbounds = [g>0; g<0];
        qred = sum(onbound(:) & exbounds(:),1)>0;
        ar.firstorderopt = norm(g(~qred));
        %         fprintf('first order optimality criterion %f (%i)\n', ar.firstorderopt, -sum(qred));
    else
        ar.firstorderopt = nan;
    end
end

if(has_error)
    rethrow(error_id)
end

if(nargout>0 && ~qglobalar)
    varargout{1} = ar;
else
    varargout = cell(0);
end

if ~silent
    arGetMerit
end

function c = my_equals(a,b)
c = a(:)==b(:);
c = c';
