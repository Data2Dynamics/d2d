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

if(nargin==0 || ~isstruct(varargin{1}))
    global ar %#ok<TLEV>
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

if(~isempty(varargin))
    sensi = varargin{1};
else
    sensi = true;
end
if(length(varargin)>1)
    fine = varargin{2};
else
    fine = false;
end
if(length(varargin)>2)
    dynamics = varargin{3}; %#ok<NASGU>
else
    dynamics = sum(ar.qDynamic == 1 & ar.qFit == 1) > 0 || ~sensi; %#ok<NASGU>
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

ar.stop = 0;

% propagate parameters
for m=1:length(ar.model)  
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
if(fine && sensi)
    ar = initFineSensis(ar);
end

% call mex function to simulate models
feval(ar.fkt, ar, fine, ar.config.useSensis && sensi, dynamics)

% integration error ?
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        if(ar.model(m).condition(c).status>0)
            error('arSimuCalc failed at %s', ar.info.arsimucalc_flags{ar.model(m).condition(c).status});
        elseif(ar.model(m).condition(c).status<0)
            error('cvodes failed at %s for model %i, condition %i', ...
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



% Initialize arrays for fine sensitivities with zeros
function ar = initFineSensis(ar)

for m = 1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d = 1:length(ar.model(m).data)
            ar.model(m).data(d).syFineSimu = zeros(length(ar.model(m).data(d).tFine), length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
            ar.model(m).data(d).systdFineSimu = zeros(length(ar.model(m).data(d).tFine), length(ar.model(m).data(d).y), length(ar.model(m).data(d).p));
        end
    end
    for c = 1:length(ar.model(m).condition)
        ar.model(m).condition(c).suFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).u), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).svFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).vs), length(ar.model(m).condition(c).p));
        ar.model(m).condition(c).sxFineSimu = zeros(length(ar.model(m).condition(c).tFine), length(ar.model(m).x), length(ar.model(m).condition(c).p));
    end
end
