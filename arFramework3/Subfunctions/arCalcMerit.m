% This function updates the objective functions used for fitting for the
% currently chosen parameters.
% 
% This function calls
%   - arSimu which in turn calls arCalcRes
%   - arCollectRes
% 
% This function replaces arChi2 and is called by arFit
% 
% Example: (call of arFit)
% ar = arCalcMerit(ar, true, ar.p(ar.qFit==1));
% 
%-----------------------------------------
%
%   Detailed description: 
%   (taken from arChi2):
% 
% arCalcMerit(sensi, pTrial, dynamics, doSimu)
%   sensi:          propagate sensitivities         [false]
%                   this argument is passed to arSimu
%   pTrial:         trial parameter of fitting
%   dynamics:       force evaluation of dynamics    [false]
%   doSimu          should arSimu be called         [true]
% 
% or
%
% ar = arCalcMerit(ar, sensi, pTrial, dynamics)
%   ar:             d2d model/data structure
% 
% Possible calls:
% arCalcMerit        then:
%               qglobalar = true
%               sensi = true
% 
% arCalcMerit(ar)    here, the global ar is overwritten by the argument
%               qglobalar = false
%               sensi = true
% 
% arCalcMerit(ar,sensi)
% arCalcMerit(sensi)
%               like arCalcMerit, but sensi can be set 
% 
% arCalcMerit(ar,sensi,ptrial)
% arCalcMerit(sensi,ptrial)
%               like arCalcMerit(ar,sensi), but 
%               silent = true
%               ar.p is set to ptrial
%               this is the only possiblity to set 'silent' to true
% 
% arCalcMerit(sensi,ptrial,dynamics)
% arCalcMerit(ar,sensi,ptrial,dynamics)
%               dynamics is passed to arSimu, otherwise arSimu is called
%               without this argument, i.e. with its default 
%               
% arCalcMerit(sensi,ptrial,dynamics,doSimu)
% arCalcMerit(ar,sensi,ptrial,dynamics,doSimu)
%           doSimu  can be set to false, then the residuals are calculated
%           without updating the model trajectories (e.g. if ar.qFit or ar.
%           model.data.qFit) has been changed.



function varargout = arCalcMerit(varargin)


global ar

nargs = nargin;
% The possiblity providing ar as an argument and to use of qglobalar==0 is
% obsolete because the gloal "ar" is overwritten anyway in arSimu 
% Implementation due to backwards compability:
if nargs>0 && isstruct(varargin{1})
    ar = varargin{1};
    varargin = varargin(2:end);
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
qPositiveX = cell(1,length(ar.model));

if(~isfield(ar.config, 'nCVRestart') || isnan(ar.config.nCVRestart))
    nCVRestart = 10;
else
    nCVRestart = ar.config.nCVRestart;
end

for i = 1:nCVRestart
    try
        if doSimu
            if(qglobalar)  % since ar is overwritten anyway in arSimu, the possiblity to use of qglobalar obsolete
                arSimu(sensi, false, dynamics);
            else
                ar = arSimu(ar, sensi, false, dynamics);
            end
        end
        has_error = false;
        break
    catch error_id
        has_error = true;
        if(~silent)
            disp(error_id.message);
        end
        if nCVRestart > 1
            if strcmp(error_id.identifier,'MATLAB:UndefinedFunction')
                break
            else
                error_printed = 0;
                for m=1:length(ar.model)
                    for c=1:length(ar.model(m).condition)
                        if(ar.model(m).condition(c).status==-1)
                            % CV_TOO_MUCH_WORK
                            if(isempty(qPositiveX{m}))
                                qPositiveX{m} = ar.model(m).qPositiveX;
                                ar.model(m).qPositiveX(:) = 0;
                            else
                                ar.config.maxsteps = (1+.2*(i-1))*maxsteps;
                                if(~error_printed)
                                    arFprintf(1, 'Integration error, restarting %d / %d with 20%% increased maxsteps.\n',i-1,nCVRestart)
                                    error_printed = 1;
                                end
                            end
                        elseif(ar.model(m).condition(c).status<-1)
                            ar.config.atol = (1+.05*i)*atol;
                            ar.config.rtol = (1+.05*i)*rtol;
                            if(~error_printed)
                                arFprintf(1, 'Integration error, restarting %d / %d with 5%% increased precision.\n',i,nCVRestart)
                                error_printed = 1;
                            end
                        else
                            error( error_id.message );
                        end
                    end
                end

            end
        end 
    end
end
ar.config.atol = atol;
ar.config.rtol = rtol;
ar.config.maxsteps = maxsteps;
for m=1:length(ar.model)
    if(~isempty(qPositiveX{m}))
        ar.model(m).qPositiveX = qPositiveX{m};
    end
end
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(isfield(ar.model(m).condition(c), 'xExpSimu'))
            if(sum((min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) ./ range(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),1) < -ar.config.rtol) & (min(ar.model(m).condition(c).xExpSimu(:,qPositiveX{m}==1),[],1) < -ar.config.atol)) > 0)
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

function y = range(x, dim)
y = max(x, [], dim) - min(x, [], dim);

function c = my_equals(a,b)
c = a(:)==b(:);
c = c';
