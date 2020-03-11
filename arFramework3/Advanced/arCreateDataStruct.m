% arCreateDataStruct(m,[pold],[fp],varargin)
% 
% This function can be used to add data structs before compiling, i.e. in
% addition to arLoadData
% 
%   pold    parameters which have to be replaced (like conditions) in order
%           to get data-specific parameters
% 
%   fp      new parameter names
% 
% 
% Example:
% D = arCreateDataStruct(1)
% 
% Example:
% D = arCreateDataStruct(1,{'maxtime'},{'100'});
% arAddDataStruct(D);
% arCompileAll;
% 
% 
% See also arCondition2NewDataStruct, arCondition2NewDataStruct2,
% arCondition2NewDataStruct3, arFitTransientFunction2

function D = arCreateDataStruct(m,pold,fp,varargin)
global ar
if ~exist('pold','var') || isempty(pold)
    pold = cell(0);
elseif ~iscell(pold)
    pold = {pold};% make cell 
end
if ~exist('fp','var') || isempty(fp)
    fp = cell(0);
elseif ~iscell(fp)
    fp = {fp}; % make cell 
end

if length(pold) ~= length(fp)
    error('length(pold) ~= length(fp)')
end

fns = varargin(1:2:(end-1));
fvals = varargin(2:2:end);

D = arEmptyDataStruct(m);

% fy, fystd, fu are inhereted from the model (like in arLoadData)
if isfield(ar.model(m),'fy')
    D.fy = ar.model(m).fy;
else
    D.fy = cell(0);
end
if isfield(ar.model(m),'fystd')
    D.fystd = ar.model(m).fystd;
else
    D.fystd = cell(0);
end
D.fu = ar.model(m).fu;


for i=1:length(fns)
    D.(fns{i}) = fvals{i};
end

%% Default values (if not yet available, i.e. if  not set via varargin)
if isempty(D.tLim)
    if ~isempty(D.tExp)
        D.tLim = [min(ar.model(m).tLim(1),min(D.tExp)),max(ar.model(m).tLim(2),max(D.tExp))];
    else
        D.tLim = ar.model(m).tLim;
    end
end
if isempty(D.tLimExp)
    D.tLimExp = [min(D.tExp),max(D.tExp)];
end
if isempty(D.logplotting)
    D.logplotting = zeros(size(D.y));
end
if isempty(D.logfitting)
    D.logfitting = zeros(size(D.y));
end
if isempty(D.normalize)
    D.normalize = zeros(size(D.y));
end
if isempty(D.tUnits) % required for plotting
    D.tUnits = {'T'  'min'  'time'};
end
if isempty(D.yUnits) % required for plotting
    D.yUnits = repmat({'C'  'au'  'conc'},size(D.y'));
end

%% Performing checks (i.e. whether varargins were reasonable
if size(D.yExp,2) ~= size(D.y,2)
    error('size(D.yExp,2) ~= size(D.y,2)');
elseif size(D.yExp,1) ~= size(D.tExp,1)
    error('size(D.yExp,2) ~= size(D.y,2)');
elseif sum(abs(size(D.yExpStd) - size(D.yExp)))~=0
    error('size(D.yExp) ~= size(D.yExpStd)');
elseif length(D.tLim) ~= 2
    error('length(D.tLim)~=2')
elseif length(D.tLimExp) ~= 2 && ~isempty(D.tLimExp)
    error('length(D.tLimExp)~=2')
elseif (sum(isnan(D.tLim)) + sum(isnan(D.tLimExp)))>0
    error('NaNs in tLim or tLimExp')
elseif sum(abs(size(D.y)-size(D.logplotting)))~=0
    error('size(D.logplotting)~=D.y');
elseif sum(abs(size(D.y)-size(D.logfitting)))~=0
    error('size(D.logfitting)~=D.y');
elseif sum(abs(size(D.y)-size(D.normalize)))~=0
    error('size(D.normalize)~=D.y');
elseif length(D.tUnits)~=3
    error('length(D.tUnits)~=3');
%elseif size(D.yUnits,2)~=3
%    error('size(D.yUnits,2)~=3');
end


%% Now collect different parameter variables:
varlist = cellfun(@symvar, D.fu, 'UniformOutput', false);
D.pu = setdiff(vertcat(varlist{:}), {ar.model(m).t, ''}); %R2013a compatible
varlist = cellfun(@symvar, D.fy, 'UniformOutput', false);
D.py = setdiff(setdiff(vertcat(varlist{:}), union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible

for j=1:length(D.fy)
    varlist = symvar(D.fy{j});
    D.py_sep(j).pars = setdiff(setdiff(varlist, union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z)), {ar.model(m).t, ''}); %R2013a compatible
    
    % exclude parameters form model definition
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).px);
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).pu);
end

varlist = cellfun(@symvar, D.fystd, 'UniformOutput', false);
D.pystd = setdiff(vertcat(varlist{:}), union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), ... %R2013a compatible
    D.y), ar.model(m).t));

for j=1:length(D.fystd)
    varlist = symvar(D.fystd{j});
	D.py_sep(j).pars = union(D.py_sep(j).pars, ... %R2013a compatible
        setdiff(varlist, union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), D.y))); %R2013a compatible
    
    % exclude parameters form model definition
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).px);
    D.py_sep(j).pars = setdiff(D.py_sep(j).pars, ar.model(m).pu);
end

varlist = cellfun(@symvar, D.fp, 'UniformOutput', false);
D.pcond = setdiff(vertcat(varlist{:}), D.p); %R2013a compatible


ptmp = union(ar.model(m).px, ar.model(m).pu);
D.p = union(ptmp, union(D.pu, D.py)); %R2013a compatible

varlist = cellfun(@symvar, D.fystd, 'UniformOutput', false);
D.pystd = setdiff(vertcat(varlist{:}), union(union(union(union(ar.model(m).x, ar.model(m).u), ar.model(m).z), D.y), ar.model(m).t));
D.pystd = setdiff(D.pystd, D.p); %Remove dynamic variables from error model parameters (like in arLoadData, it rarely occurs but might be a problem for fitting)
D.p = union(D.p, D.pystd); %R2013a compatible
if ( size( D.p, 1 ) ~= 1 )
    D.p = D.p.';
end

D.fp = transpose(D.p);

% execute substitutions from ar.model.p/ar.model.fp
[int, iA, iB] = intersect(D.fp, ar.model.p);
D.fp(iA) = ar.model.p(iB);

%% now replace fp which was provided as function argument:
[~,ia,ib] = intersect(D.p,pold);
for i=1:length(ia)
    D.fp{ia(i)} = fp{ib(i)};
end
ia = find(contains(D.fy,pold));
for i=1:length(ia)
    for j=1:length(pold)
        D.fy{i} = strrep(D.fy{ia(i)},pold{j},fp{j});
    end
end
for i=1:length(D.fystd)
    for j=1:length(pold)
        D.fystd{i} = strrep(D.fystd{i},pold{j},fp{j});
    end
end


%% This function creates an empty data struct
function D = arEmptyDataStruct(m)
global ar
if ~exist('m','var') 
    m = [];
end

if isempty(m) || ~isfield(ar.model(m),'data') || isempty(ar.model(m).data)
                       D = struct;
             D.condition = [];
           D.description = {'Created with arEmptyDataStruct.m'};
          D.doseresponse = false;
                    D.fp = cell(0);
                    D.fu = cell(0);
                    D.fy = cell(0);
                 D.fystd = cell(0);
            D.logfitting = [];
           D.logplotting = [];
                  D.name = 'arEmptyDataStruct';
             D.normalize = [];
                     D.p = cell(0);
                  D.path = pwd;
                 D.pcond = cell(0);
                 D.prand = cell(0);
                    D.pu = cell(0);
                    D.py = cell(0);
                D.py_sep = struct;
                 D.pystd = cell(0);
             D.rand_type = [];
    D.response_parameter = '';
                     D.t = 'time';
                  D.tExp = [];
                  D.tLim = [];
               D.tLimExp = [];
                D.tUnits = cell(0);
                D.uNames = cell(0);
                     D.y = cell(0);
                  D.yExp = [];
               D.yExpRaw = [];
               D.yExpStd = [];
            D.yExpStdRaw = [];
                D.yNames = cell(0);
                D.yUnits = cell(0);
%               D.checkstr = '';
%                    D.fkt = '';
%                  D.cLink = [];
                 
else % use model(m).data(1) as template:
    D = ar.model(m).data(1); 
    fn = fieldnames(D);
    for f=1:length(fn)
        if iscell(D.(fn{f}))
            D.(fn{f}) = cell(0);
        elseif isnumeric(D.(fn{f}))
            D.(fn{f}) = [];
        elseif islogical(D.(fn{f}))
            D.(fn{f}) = [];
        elseif ischar(D.(fn{f}))
            D.(fn{f}) = '';
        elseif isstruct(D.(fn{f}))
            D.(fn{f}) = struct;
        else
            error('This case is not yet implemented.');
        end
    end
    
end