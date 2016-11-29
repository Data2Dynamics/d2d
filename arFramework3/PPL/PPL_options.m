function options = PPL_options(varargin)
global ar
% Set options for calculation of prediction profile likelihood

% gammas : vector of correction strengths for computation of integration
% steps for the vector of states ix (default is 0.1 for all ix)
%
% stepsize : Integration steps in time, e.g. 1, 0.5 ... (default is tFine stepsize)
%
% onlyProfile : Don't integrate prediction/validation confidence line over
% time, just compute profiles for specified t (default is false)
%
% tEnd : optional stopping time for the confidence band integration
% (default is tLim(end))
%
% rel_increase : increase of xTrial in computation of full profile for t's
% (default is 0.15)
%
% ed_steps : take integration steps or go parallel to orginial state solution (default is true)
%
% fineInt : do precise integration with a correction cycle at every
% integration step (default is true)
%
% backward : perform integration backwards in time (default is false)
%
% xstd : standard deviation of additional data point. If profile on data
% and ar.config.fiterrors=1, the fitted errors of the respective data is
% used (default is 0.1)
%
% doPPL : integrate prediction bands without assumed error on measurements
% (? confidence bands and is based on prediction profiles of Clemens paper
% (default is false ? validation profiles))
%
% n_start : total steps in computation of full profile per side (default is 100)
%
% whichT : which entry of the provided time-vector is taken as starting
% point for integration
%
% dir : 1 for only upper CI profile, -1 for lower, 0 for both (default is 0)
%
% alpha_level : Confidence level for profile (default 0.05)
%
% Helge Hass, 2014 (helge.hass@fdm.uni-freiburg.de)

if (nargin == 0) && (nargout == 0)
  fprintf('          Gammas: [ vector of correction strengths {1/stepsize} ]\n');
  fprintf('        stepsize: [ integration step size {tFine} ]\n');
  fprintf('     onlyProfile: [ true | {false} ]\n');
  fprintf('            tEnd: [ Alternative integration endpoint ]\n'); 
  fprintf('    rel_increase: [ %% of x increase in profile calculation {0.15} ]\n');
  fprintf('        ed_steps: [ {true} | false ]\n');
  fprintf('         fineInt: [ {true} | false ]\n');  
  fprintf('        backward: [ true | {false} ]\n');
  fprintf('            xstd: [ standard dev of auxiliary data point {0.1} ]\n');
  fprintf('           doPPL: [ true = prediction or false = validation profiles {false} ]\n');
  fprintf('         n_start: [ steps of profile {100} ]\n');
  fprintf('          whichT: [ seed for integration (in time vector) {1} ]\n');
  fprintf('             dir: [ only upper/lower profile? {0} ]\n');
  fprintf('     alpha_level: [ Confidence level {0.05} ]\n');
  fprintf('\n');
  return;
end

Names = [
    'gammas          '
    'stepsize        '
    'onlyProfile     '
    'tEnd            '
    'rel_increase    '
    'ed_steps        '             
    'fineInt         '
    'backward        '
    'xstd            '
    'doPPL           '
    'n_start         ' 
    'whichT          '
    'dir             '
    'alpha_level     '
    ];

%Set some default options
if(~isfield(ar.ppl,'options'))
    ar.ppl.options.onlyProfile = false;
    ar.ppl.options.rel_increase = 0.15;
    ar.ppl.options.ed_steps = true;
    ar.ppl.options.fineInt = true;
    ar.ppl.options.backward = false;
    ar.ppl.options.xstd = 0.1;
    ar.ppl.options.doPPL = false;
    ar.ppl.options.n_start = 100;
    ar.ppl.options.whichT = 1;
    ar.ppl.options.dir = 0;
    ar.ppl.options.alpha_level = 0.05;
end

m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
% options = [];
% for j = 1:m
%   options.(deblank(Names(j,:))) = [];
% end
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(message('No proper name %s', arg));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        ar.ppl.options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error(message('No argument for a given option'));
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~ischar(arg)
      error(message('Option must be a char', arg));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(message('Could not find option %s', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
            matches = deblank(Names(j(1),:));
        for k = j(2:length(j))'
                matches = [matches ', ' deblank(Names(k,:))]; %#ok<AGROW>
        end
            error(message('Found more than one match for option %s',arg));
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    ar.ppl.options.(deblank(Names(j,:))) = arg;
    if(strcmp(deblank(Names(j,:)),'xstd'))
        ar.ppl.xstd_auto = 1;
    end
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(message('No value found for option %s', arg));
end

options = 'Did set PPL options! \n';
